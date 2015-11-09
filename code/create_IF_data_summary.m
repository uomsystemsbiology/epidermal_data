%% create_IF_data_summary.m
% This MATLAB script reproduces the summary figures shown in Additional
%  File 4 of the GigaScience Data Note:
%   Cursons et al. (2015). Spatially-transformed fluorescence image data 
%    for ERK-MAPK and selected proteins within human epidermis.
%    GigaScience. Submitted Sept 2015.
%   http://dx.doi.org/doi-not-yet-known
%
% It is designed to be executed upon the corresponding epidermal 
%  fluorescence data published with this Data Note:
%   Cursons, J; Angel, C, E; Hurley, D, G; Print, C, G; Dunbar, P; Jacobs, 
%    M, D; Crampin, E, J (2015): Supporting data for "Spatially-transformed
%    fluorescence image data for ERK-MAPK and selected proteins within 
%    human epidermis". GigaScience Database.
%   http://dx.doi.org/10.5524/100168 
%
% A number of functions are used by this script, some of which have
%  dependencies upon MATLAB Toolboxes:
%   - Image Processing Toolbox: this script directly calls the imread, 
%           tforminv, and cp2tform functions
%   - Neural Network Toolbox: ?
%   - Statistics and Machine Learning Toolbox: ?
%
% This script also uses the 'rotate_image' function from the MATLAB File
%   Exchange: File ID: #4071
%       http://www.mathworks.com/matlabcentral/fileexchange/4071-rotate-image
%       rotate_image.m was created by Ohad Gal: 
%           http://www.mathworks.com/matlabcentral/profile/authors/869576-ohad-gal
%
% Please note that this script is a behemoth and could be (may be in the
%  future) written using a number of functions; however there are a number
%  of complexities associated with passing figure handles in/out of
%  functions.
% Also note that these scripts were developed when I was using the data for
%  network inference methods so I regularly refer to a 'node' which is
%  shorthand for a protein within a specific-subcellular localisation
%
% This script was created by Joe Cursons at the University of Melbourne
%   Systems Biology Laboratory:
%       joseph.cursons@unimelb.edu.au
%
% Last Updated: 09/11/15
%
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Input Parameters - User Defined Settings
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %   
%in this section I specify a number of arrays which can be used for 
% subsequent indexing between different data arrays to create the output
% figures

%examine loess curves/smoothed outputs for all three patients contained in
% these data 
arrayPatients = 1:3;
numPatients = size(arrayPatients,2); 
 
%as noted above, these scripts were developed when I was using network
% methods and the term node refers to a sub-cellular localised protein; in
% total there are 17 target proteins over three possible sub-cellular 
% localisations (some are "empty"/no observation) --> 51 total nodes
numNodesTotal = 51; 

%arrayGroupedNodes specifies sub-cellular localisations as columns and
% proteins as nodes, containing index numbers for the 'nodes' (this array
% is in-part a relic of a previous script structure)
arrayGroupedNodes = [  1,  0, 35;                   % ITGB1
                       3,  0, 37;                   % ITGB4
                       4,  0,  0;                   % SFN
                       5, 22, 39;                   % CALM
                       6, 23,  0;                   % Raf
                       7, 24,  0;                   % pRaf
                       8, 25,  0;                   % MEK
                       9, 26,  0;                   % pMEK
                       10, 27,  0;                  % ERK
                       11, 28,  0;                  % pERK
                       12, 29,  0;                  % JunB
                       13, 30,  0;                  % cJun
                       14, 31,  0;                  % cFos
                       15, 32,  0;                  % Fra2
                       16,  0,  0;                  % K10
                       17,  0,  0  ];               % K14
numOutputProteins = size(arrayGroupedNodes,1);

%arrayOutputNodeOrder provides a condensed version of arrayGroupedNodes
% witout padding zeroes
arrayOutputNodeOrder = [ 1, 3:17, 22:32, 35, 37, 39 ];
numOutputNodes = length(arrayOutputNodeOrder);

%arrayOutputProtFolders contains protein folder names in an order which
% corresponds to the rows of arrayGrouped Nodes
arrayOutputProtFolders = { 'ITGB1'; 'ITGB4'; 'SFN'; 'CALM'; 'Raf'; ...
                           'Raf_ph'; 'MEK'; 'MEK_ph'; 'ERK'; 'ERK_ph'; ...
                           'JUNB'; 'JUNC'; 'FOSC'; 'FRA2'; 'K10'; ...
                           'K14' };

%arrayOutputProtFolders contains LaTeX compatible strings ('$$' around
% special characters) for overlaying on the figures.
arrayProteinLateXString = { ['$${\beta}1$$-Integrin $$\;$$']; ...
                            ['$${\beta}4$$-Integrin $$\;$$']; ...
                            ['14-3-3$${\sigma}$$']; ...
                            ['Calmodulin $$\;$$']; ...
                            ['Raf-1']; ...
                            ['phospho-Raf-1 $$\; \; \; \;$$']; ...
                            ['MEK-1/2 $$\;$$']; ...
                            ['phospho-MEK-1/2 $$\; \; \; \;$$']; ...
                            ['ERK-1/2 $$\;$$']; ...
                            ['phospho-ERK-1/2 $$\; \; \; \;$$']; ...
                            ['Jun-B']; ...
                            ['c-Jun']; ...
                            ['c-Fos']; ...
                            ['Fra-2']; ...
                            ['Keratin-10 $$\;$$']; ...
                            ['Keratin-14 $$\;$$'] };
                        
                        
%specify the patient for the 'representative image data' of the protein
arrayPatientToDisplay = [ 1; 2; 2; 1; 2; ...
                          2; 1; 2; 2; 1; 3; ...
                          2; 3; 2; 3; 1 ];
%specify the z-position/image to use from the patient specified above                      
arrayRepImageZ = { 63; 15; 21; 25; 15; ...
                    5; 14; 5; 30; 6; 20; ...
                    10; 2; 5; 25; 25 };

%the color scheme Red = Pat1, Green = Pat2 and Blue = Pat3 is used across
% these data
arrayPatientColors = { [ 1.0 0.0 0.0 ];
                       [ 0.0 1.0 0.0 ];
                       [ 0.0 0.0 1.0 ] };
%use 'muted' R/G/B colors for the data cloud 
arrayDataCloudColor = [ 1.0, 0.8, 0.8;
                        0.8, 1.0, 0.8;
                        0.8, 0.8, 1.0 ];
                   

%specify the spatial discretisation settings for the epidermal tissue
numSpatialBins = 28;
arrayNormDistRatio = [ 1 4 2 ];
numTissueLayers = size(arrayNormDistRatio, 2);

%specify the smoothing parameter for the loess smoothing (the relative
% fraction of 'nearest' data used for averaging)
numLoessWindowSize = 0.5;

%specify the size of the confidence interval for the loess smoothing
numUniDirConfInt = 0.90;
                   
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Input Parameters - Figure Display Settings
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %   
%set the scaling factor for displaying the image on the screen
numFigScaleMult = 4;

%specify font and marker sizes etc within the output images
numPlotFontSizeTitle = 11;
numPlotFontSize = 9;
numPlotDataCloudMarkerSize = 2;
numPlotLowessMarkerSize = 8;
numOverlayBorderWidth = 2;
numSubFigLabelFontSize = 18;

%specify the size of the input images
numImageSizePixels = 1024;

%specify the number of pixels within 10 um for each of the image data sets
% for drawing the scale bar
arrayNumPixelsIn10um = [ 067; 065; 069; 074; 058; ...
                         073; 072; 058; 063; 069; 077; ...
                         065; 090; 061; 079; 067 ];
                     
%specify the scale bar relative width and spacing values (from the image
% border)
numFractionScaleBarSpacing = 0.1;
numFractScaleBarWidth = 0.025;

%specify the number of images in each image stack
%TODO: FINISH POPULATING THIS
arrayImagesInStack = [ 66; 64; 32; 64; 00; ...
                     00; 00; 00; 00; 00; 00; ...
                     00; 00; 00; 71; 00; ];
                 
%specify the angle for rotating the display/representative image
arrayRepImageRotate = [ 180; 180; 090; 180; 000; ...
                        180; 180; 000; 000; 180; 180; ...
                        000; 000; 180; 090; 180 ];
                    
%specify the angle for rotating the zoomed display/representative image
arrayZoomedRepImageRotate = [ 000; -20; -15; -15; -20; ...
                              -10; 5; -20; 000; -20; -05; ...
                              -20; -45; -30; 030; 000  ];
                          
%specify the x-y coordinates for the zoomed region (before rotation)
arrayZoomedRepImageCrop = [  311   625   296   855;        % ymin ymax xmin xmax
                             641  1045   276   995;
                             201   605   246   965;
                             566   880   466  1025;
                             291   695   346  1065;
                             481   840   236   975;
                             321   635   361   920;
                             351   755   276   995;
                             246   750    66   945;
                             431   745   526  1085;
                             631   900   306   785;
                             401   805   356  1075;
                             461   820   416  1055;
                             411   815   316  1035;
                             456   680   516   915;
                             301   705   216   935 ];

%specify the 'figure type' for the different images - there are different
% numbers of sub-cellular localisations, so different plot types are
% included etc, this is achieved by specifying these values over a number
% of arrays
%the 'other image' types for display (excluding the representative image)
% 2 = surface plot, 3 = 3D render
arrayOtherImageTypes = cell(numPatients,1); 
%the angle for rotation of this other image
arrayOtherImageRot = zeros(numPatients,2,'double');
%the x-y co-ordinates for the 'other image'
% xmin xmax ymin ymax zmin zmax: 0000 if not defined
arrayOtherImageCoOrds = cell(numPatients,2); 
%the 'view' parameter ("camera angle") for the other image
arrayOtherImageView = cell(numPatients,2);
%the number of times to perform interpolation for surface plots
arrayNumSurfaceInterp = ones(numPatients,2,'uint8');
%the number of isosurface groups to use for volume/3D render 'other images'
arrayIsoSurfGroups = cell(numPatients,1);

%ITGB1 - surface plot and a 3D render
arrayOtherImageTypes{1} = [ 2 3 ];
arrayOtherImageCoOrds{1,1} = [ 0334 0525 0495 0605 0063 0063 ];
arrayOtherImageCoOrds{1,2} = [ 0334 0525 0495 0605 0045 0065 ];
arrayOtherImageRot(1,:) = [ 0, 0 ];
arrayNumSurfaceInterp(1,:) = [ 4, 1 ];
arrayOtherImageView{1,1} = [ -18 80 ];
arrayOtherImageView{1,2} = [ -30 44 ];
arrayIsoSurfGroups{1} = {[ 099 129 ];
                         [ 130 169 ];
                         [ 170 255 ]};

%ITGB4 - surface plot and a 3D render
arrayOtherImageTypes{2} = [ 2 3 ];
arrayOtherImageCoOrds{2,1} = [ 0325 0460 0300 0400 0015 0015 ];
arrayOtherImageCoOrds{2,2} = [ 0325 0460 0300 0400 0015 0035 ];
arrayOtherImageRot(2,:) = [ -20, -20 ];
arrayOtherImageView{2,1} = [ -19 58 ];
arrayOtherImageView{2,2} = [ -19 58 ];
arrayIsoSurfGroups{2} = {[ 050 069 ];
                         [ 070 099 ];
                         [ 100 255 ]};

%SFN - surface plot
arrayOtherImageTypes{3} = [ 2 ];
arrayOtherImageCoOrds{3,1} = [ 0330 0520 0300 0465 0000 0000 ];
arrayOtherImageRot(3,:) = [ -20, 0 ];
arrayNumSurfaceInterp(3,:) = [ 4, 1 ];
arrayOtherImageView{3,1} = [ -18 78 ];               
                    
%CALM - surface plot and a 3D render
arrayOtherImageTypes{4} = [ 2 3 ];
arrayOtherImageCoOrds{4,1} = [ 0495 0685 0740 0865 0025 0025 ];
arrayOtherImageCoOrds{4,2} = [ 0495 0685 0740 0865 0015 0035 ];
arrayOtherImageRot(4,:) = [ -15, -15 ];
arrayNumSurfaceInterp(4,:) = [ 4, 1 ];
arrayOtherImageView{4,1} = [ -17 72 ];
arrayOtherImageView{4,2} = [ -20 60 ];
arrayIsoSurfGroups{4} = {[ 075 099 ];
                         [ 100 119 ];
                         [ 120 255 ]};
                    
%Raf - surface plot
arrayOtherImageTypes{5} = [ 2 ];
arrayOtherImageCoOrds{5,1} = [ 0400 0520 0550 0615 0000 0000 ];
arrayOtherImageRot(5,:) = [ 0, 0 ];
arrayNumSurfaceInterp(5,:) = [ 4, 1 ];
arrayOtherImageView{5,1} = [ -21 88 ];

%Raf_ph - surface plot
arrayOtherImageTypes{6} = [ 2 ];
arrayOtherImageCoOrds{6,1} = [ 0450 0620 0510 0650 0000 0000 ];
arrayOtherImageRot(6,:) = [ 0, 0 ];
arrayNumSurfaceInterp(6,:) = [ 4, 1 ];
arrayOtherImageView{6,1} = [ -21 86 ];

%MEK - surface plot
arrayOtherImageTypes{7} = [ 2 ];
arrayOtherImageCoOrds{7,1} = [ 0565 0670 0635 0760 0000 0000 ];
arrayOtherImageRot(7,:) = [ 0, 0 ];
arrayNumSurfaceInterp(7,:) = [ 4, 1 ];
arrayOtherImageView{7,1} = [ -21 86 ];

%MEK_ph - surface plot
arrayOtherImageTypes{8} = [ 2 ];
arrayOtherImageCoOrds{8,1} = [ 0445 0670 0595 0750 0000 0000 ];
arrayOtherImageRot(8,:) = [ -20, 0 ];
arrayNumSurfaceInterp(8,:) = [ 4, 1 ];
arrayOtherImageView{8,1} = [ -21 86 ];

%ERK - surface plot
arrayOtherImageTypes{9} = [ 2 ];
arrayOtherImageCoOrds{9,1} = [ 0345 0470 0430 0575 0000 0000 ];
arrayOtherImageRot(9,:) = [ 0, 0 ];
arrayNumSurfaceInterp(9,:) = [ 4, 1 ];
arrayOtherImageView{9,1} = [ -16 82 ];

%ERK_ph - surface plot
arrayOtherImageTypes{10} = [ 2 ];
arrayOtherImageCoOrds{10,1} = [ 0660 0780 0440 0540 0000 0000 ];
arrayOtherImageRot(10,:) = [ 0, 0 ];
arrayNumSurfaceInterp(10,:) = [ 4, 1 ];
arrayOtherImageView{10,1} = [ -15 84 ];

%JUNB - surface plot
arrayOtherImageTypes{11} = [ 2 ];
arrayOtherImageCoOrds{11,1} = [ 0460 0580 0750 0865 0000 0000 ];
arrayOtherImageRot(11,:) = [ -5, 0 ];
arrayNumSurfaceInterp(11,:) = [ 4, 1 ];
arrayOtherImageView{11,1} = [ -16 86 ];

%JUNC - surface plot
arrayOtherImageTypes{12} = [ 2 ];
arrayOtherImageCoOrds{12,1} = [ 0565 0655 0595 0690 0000 0000 ];
arrayOtherImageRot(12,:) = [ -20, 0 ];
arrayNumSurfaceInterp(12,:) = [ 4, 1 ];
arrayOtherImageView{12,1} = [ -16 86 ];

%FOSC - surface plot
arrayOtherImageTypes{13} = [ 2 ];
arrayOtherImageCoOrds{13,1} = [ 0530 0710 0580 0700 0000 0000 ];
arrayOtherImageRot(13,:) = [ -45, 0 ];
arrayNumSurfaceInterp(13,:) = [ 4, 1 ];
arrayOtherImageView{13,1} = [ -18 82 ];

%FRA2 - surface plot
arrayOtherImageTypes{14} = [ 2 ];
arrayOtherImageCoOrds{14,1} = [ 0650 0865 0660 0865 0000 0000 ];
arrayOtherImageRot(14,:) = [ -45, 0 ];
arrayNumSurfaceInterp(14,:) = [ 4, 1 ];
arrayOtherImageView{14,1} = [ -18 86 ];

%K10 - surface plot
arrayOtherImageTypes{15} = [ 2 ];
arrayOtherImageCoOrds{15,1} = [ 0635 0730 0400 0520 0000 0000 ];
arrayOtherImageRot(15,:) = [ 30, 0 ];
arrayNumSurfaceInterp(15,:) = [ 1, 1 ];
arrayOtherImageView{15,1} = [ 16 64 ];

%K14 - surface plot
arrayOtherImageTypes{16} = [ 2 ];
arrayOtherImageCoOrds{16,1} = [ 0245 0490 0200 0330 0000 0000 ];
arrayOtherImageRot(16,:) = [ 0, 0 ];
arrayNumSurfaceInterp(16,:) = [ 4, 1 ];
arrayOtherImageView{16,1} = [ -25 76 ];


 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Path Manipulations
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  

%depending upon the OS, use a different folder separater (forward vs back
% slash)
if ispc,
    strFolderSep = '\';
elseif isunix,
    strFolderSep = '/';
else
    disp(['warning: cannot determine the operating system, defaulting to ' ...
            'forward slash for the folder path separator' ] );
    strFolderSep = '/';
end

%extract the current directory
strCurrDir = cd;

%the functions used by this script should be located in the /functions/
% folder, so add this to the MATLAB file path
if ~isdeployed
    addpath(genpath(strCurrDir));
else
    %if the machine is deployed, include the ctfroot (must be defined when
    % running the matlab compiler) in the file path
    addpath(genpath([ctfroot strFolderSep 'code']))
end


%manipulate the file path to determine the appropriate folders
if arrayCurrDirFoldSepPos(end) == length(strCurrDir),
    %there is a backslash at the end
    strBaseDir = strCurrDir(1:(arrayCurrDirFoldSepPos(end-1)));
else
    strBaseDir = strCurrDir(1:(arrayCurrDirFoldSepPos(end)));
end

%determine the relative image data folder path
strImgDataFolder = [ strBaseDir 'image' strFolderSep ];

%determine the relative processed data folder path
strProcDataFolder = [ strBaseDir 'processed' strFolderSep ];
    
%for now just output to the base directory
stringOutputDataFolder = strBaseDir;

 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Perform Pre-Processing - Extract Confidence Intervals for Residuals
%   From the Loess Curve 
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  

%extract the full spatial discretisation arrays
[arrayIntervalCentres, arrayDistIntervals, arrayDivisionIndices] = calculateSpatialDivisions( numSpatialBins, arrayNormDistRatio );
arrayLayerBoundaries = arrayDivisionIndices-1;

%calculate the confidence interval bounds for loess smoothing at various
% positions along the normalised distance coordinate
arrayLoessCIBounds = cell(numPatients, numNodesTotal, numTissueLayers);
arrayLoessCIXPos = cell(numPatients, numNodesTotal, numTissueLayers);

%move through all specified output nodes to calculate the confidence
% intervals around the loess-smoothed curve
disp(['calculating confidence intervals for residuals from the Loess curve..']);
for iNode = 1:numOutputNodes,
    
    disp([char(9) '.. for localised target ' num2str(iNode) ' of ' num2str(numOutputNodes) ]);
    
    %determine the target protein/localisation
    numNode = arrayOutputNodeOrder(iNode);
    [numOutputProt, numLoc] = find(arrayGroupedNodes == numNode);
    
    %use this to specify the data sub-folder
    if numLoc == 1,
        strLoc = 'C';
    elseif  numLoc == 2,
        strLoc = 'N';
    elseif  numLoc == 3,
        strLoc = 'M';
    else
        disp(['warning: the sample location cannot be determined from "arrayGroupedNodes" - numLoc == ' num2str(numLoc)]);
    end
    
    %step through every patient
    for iPatient = 1:numPatients,
        
        %determine the data sub-folder from the target/patient
        stringSpecificProcDataFolder = [strProcDataFolder arrayOutputProtFolders{numOutputProt} strFolderSep 'Pat_' num2str(iPatient) strFolderSep strLoc strFolderSep ];
        
        %load the analysed sample data which have been converted to
        % normalised distance
        structSampleAnalysis = loadSampAnalysis( stringSpecificProcDataFolder, '.mat' );

        %extract the signal intensity data into arrays that can be used for
        % calculating the confidence intervals around the loess curves
        %NB: to ensure backwards compatibility with older processed data
        % sets, this extracts 'objects' (grouped samples) and then attempts
        % to unpack them; in the new/latest processed data files each
        % object only contains one sample
        numObjects = length(structSampleAnalysis);
        numSamplesPerObject = size(structSampleAnalysis(1).NormDist,1);
        numPixelsPerSample = size(structSampleAnalysis(1).SigInt,2);
        %dependent upon the extracted feature information, create a
        % corresponding output array for sample positions and signal
        % intensity values
        arraySampleValues = zeros(numObjects*numSamplesPerObject, numPixelsPerSample, 'uint8');
        arraySamplePositions = zeros(numObjects*numSamplesPerObject,1, 'double');
        %populate these arrays
        for iObject = 1:numObjects,
            for iSample = 1:numSamplesPerObject,
                arraySamplePositions(  (iObject-1)*numSamplesPerObject + iSample  ) = structSampleAnalysis(iObject).NormDist(iSample);
                arraySampleValues(  (iObject-1)*numSamplesPerObject + iSample, :  ) = structSampleAnalysis(iObject).SigInt(iSample,:);
            end 
        end
        
        %rescale these data along the normalised distance axis
        arraySamplePositions = rescaleNormDist(arraySamplePositions, numSpatialBins, arrayNormDistRatio);
        
        %load the loess-smoothed data
        arrayLoessData = loadLoessCurve(stringSpecificProcDataFolder, '.mat');
        
        %extract into arrays for manipulating
        arrayTempLoessXPos = arrayLoessData(1,:);
        arrayTempLoess = arrayLoessData(2,:);
        
        %move through each tissue layer (the loess smoothing is
        % discontinuous over tissue layer boundaries/performed within
        % tissue layers)
        for iTissueLayer = 1:numTissueLayers,
            
            %identify the minimum and maximum x-positions within the tissue
            % layer
            numStartX = double(arrayDivisionIndices(iTissueLayer)-1);
            numEndX = double(arrayDivisionIndices(iTissueLayer+1)-1);
            
            %identify loess-smoothed and sampled data points within these
            % bounds
            if iTissueLayer < numTissueLayers,
                arraySamplesInTLIndex = find( (arraySamplePositions >= numStartX) & (arraySamplePositions < numEndX) );
                arrayLoessInTLIndex = find( (arrayTempLoessXPos >= numStartX) & (arrayTempLoessXPos < numEndX) );
            elseif iTissueLayer == numTissueLayers,
                arraySamplesInTLIndex = find( (arraySamplePositions >= numStartX) & (arraySamplePositions <= numEndX) );
                arrayLoessInTLIndex = find( (arrayTempLoessXPos >= numStartX) & (arrayTempLoessXPos <= numEndX) );
            end
            
            %extract just the data points within the tissue layer, and 
            % ensure that they are contiguous along the normalised distance
            % co-ordinate
            arraySamplePositionsInTL = arraySamplePositions(arraySamplesInTLIndex);
            [ arraySortedXPos, arrayIndexSortedXPos ] = sort(arraySamplePositionsInTL);
                        
            %extract just the loess-curve points within the tissue layer
            arrayLoessinTL = arrayTempLoess(arrayLoessInTLIndex);
            arrayLoessinTLXPos = arrayTempLoessXPos(arrayLoessInTLIndex);
            
            %specify the size of the array within the "output" cell array
            arrayLoessCIBounds{iPatient, numNode, iTissueLayer} = zeros(length(arrayLoessinTLXPos),2,'double');
            arrayLoessCIXPos{iPatient, numNode, iTissueLayer} = zeros(length(arrayLoessinTLXPos),1,'double');
        
            numClosestPoints = int32(double(length(arraySamplesInTLIndex))*numLoessWindowSize);
            
            %move through each individual x-point
            for iXPos = 1:length(arrayLoessInTLIndex),

                %extract the spatial position
                numXPos = arrayLoessinTLXPos(iXPos);
                arrayLoessCIXPos{iPatient, numNode, iTissueLayer}(iXPos) = numXPos;
                
                %identify index values
                numLoessVal = arrayLoessinTL(iXPos);

                %find the points closest to this x-position
                arrayXPosDiff = arraySortedXPos - numXPos;
                arrayXPosAbsDiff = abs(arrayXPosDiff);
                [ arraySortedXPosAbsDiff, arrayIndexSortedXPosAbsDiff ] = sort(arrayXPosAbsDiff);
                
                %create a pointer array for the 'closest points'
                arrayTempPointer = arrayIndexSortedXPos(arrayIndexSortedXPosAbsDiff(1:numClosestPoints));
                arrayTruePointer = arraySamplesInTLIndex(arrayTempPointer);
                
                %extract just these 'closest points' and calculate the
                % residual from the loess point
                arrayLocalSampleValues = arraySampleValues(arrayTruePointer,:);
                arrayLocalSampleDiff = double(arrayLocalSampleValues(:)) - numLoessVal;
                
                %if there are data points above/below the loess curve
                arrayAboveLoessIndex = find(arrayLocalSampleDiff > 0);
                arrayBelowLoessIndex = find(arrayLocalSampleDiff < 0);
                
                %determine the confidence interval "above" the loess curve
                if ~isempty(arrayAboveLoessIndex),
                    [arrayTempFreq, arrayTempInt] = ecdf(arrayLocalSampleDiff(arrayAboveLoessIndex));
                    arrayLoessCIBounds{iPatient, numNode, iTissueLayer}(iXPos,1) = interp1(arrayTempFreq,arrayTempInt, numUniDirConfInt);
                else
                    arrayLoessCIBounds{iPatient, numNode, iTissueLayer}(iXPos,1) = NaN;
                end
                
                %determine the confidence interval "below" the loess curve
                if ~isempty(arrayBelowLoessIndex),
                    [arrayTempFreq, arrayTempInt] = ecdf(arrayLocalSampleDiff(arrayBelowLoessIndex));
                    arrayLoessCIBounds{iPatient, numNode, iTissueLayer}(iXPos,2) = interp1(arrayTempFreq,arrayTempInt, (1-numUniDirConfInt));
                else
                    arrayLoessCIBounds{iPatient, numNode, iTissueLayer}(iXPos,2) = NaN;
                end
                
            end

        end
 
    end
    
    disp([char(9) char(9) 'approx ' num2str((iNode/numOutputNodes)*100, '%03.1f') '% complete']);
    
end

 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Create the Output Figures
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  

%move through all of the output proteins - one figure per target
disp(['creating output figures..']);
for iOutputProtein = 1:numOutputProteins,
    
    disp([ char(9) '.. for ' arrayOutputProtFolders{iOutputProtein} '(Protein ' num2str(iOutputProtein) ' of ' num2str(numOutputProteins) ')' ]);
    
    %determine output parameters for this protein
    numLocalisationsForProtein = length(find(arrayGroupedNodes(iOutputProtein,:)));
   
    %determine the folder path for the cytoplasmic sampled data 
    stringRepImageDataFolder = [ strProcDataFolder arrayOutputProtFolders{iOutputProtein} strFolderSep 'Pat_' num2str(arrayPatientToDisplay(iOutputProtein)) strFolderSep ];
    stringRepImageCytoDataFolder = [ stringRepImageDataFolder 'C' strFolderSep ];
    
    %extract the cytoplasmic sampled data    
    structSampleLocData = loadSampLocs(stringRepImageCytoDataFolder, '.mat');

    %examine the sample data locations to determine the z-positions at 
    % which samples have been collected
    arrayAllSliceZPositions = zeros(length(structSampleLocData),1,'uint8');
    for iSampleGroup = 1:length(structSampleLocData),
        arrayAllSliceZPositions(iSampleGroup) = structSampleLocData(iSampleGroup).ZPosition;
    end
    %make this list unique
    arrayZSlicesSampled = unique(arrayAllSliceZPositions);
        
    %determine the folder path for the representative image
    stringSpecificImageDataFolder = [strImgDataFolder arrayOutputProtFolders{iOutputProtein} strFolderSep 'Pat_' num2str(arrayPatientToDisplay(iOutputProtein)) strFolderSep 'image_data' strFolderSep ];   
    
    %examine contents of the folder containing the representative image
    arrayDirContents = dir(stringSpecificImageDataFolder);
    %and identify those files containing the string "_Series"; identifying
    % them as components of a fluorescence z-stack
    arraySeriesImageFlag = false(length(arrayDirContents), 1);
    for iFile = 1:length(arrayDirContents),
        if ~isempty(strfind(arrayDirContents(iFile).name, '_Series')),
            arraySeriesImageFlag(iFile) = true;
        end
    end
    %identify all images belonging to this series/stack
    numFirstSeriesImageIndex = find(arraySeriesImageFlag, 1, 'first');
    %and identify the full file name for the first image in the stack
    stringImageSeriesName = arrayDirContents(numFirstSeriesImageIndex).name;
    
    %determine the file path for the series image
    stringInputImagePath = [stringSpecificImageDataFolder stringImageSeriesName];
    
    %examine the structure of the string specifying the series image, and
    % determine the characters that correspond to the z-position
    numImagePathZPosition = strfind(stringInputImagePath, '_z0');
    
    %load the image
    stringOutputZ = num2str(arrayRepImageZ{iOutputProtein}, '%02u');
    stringInputImagePath(numImagePathZPosition+3:numImagePathZPosition+4) = stringOutputZ;
    imageRep = imread(stringInputImagePath);
       
    
    %rotate the representative image by the angle specified at the
    % beginning of this script
    %NB: this takes the input image as an array of doubles
    [ imageRepRot, ~ ] = rotate_image( arrayRepImageRotate(iOutputProtein), double(imageRep), [1; 1] );
    
    %convert back to an 8-bit unsigned integer
    imageRepRot = uint8(imageRepRot);
    
    %rescale the brightness to improve printed display and save as a new
    % image file (i.e. don't modify the original image data)
    numMinBrightnessThresh = 5;
    numMaxBrightnessScaleCoefficient = 1;
    numNonLinearTransCoefficient = 0.8;
    imageRepBright = alterPixelIntensityForDisplay( imageRepRot, numMinBrightnessThresh, numNonLinearTransCoefficient, numMaxBrightnessScaleCoefficient );

    %convert to RGB
    imageRepOut = cat(3, imageRepBright*0, imageRepBright, imageRepBright*0);
    
    %draw a scale bar at the bottom right
    numScaleBarHorStop = uint16(size(imageRepBright,2)*(1-numFractionScaleBarSpacing));
    numScaleBarHorStart = uint16(numScaleBarHorStop - arrayNumPixelsIn10um(iOutputProtein));
    numScaleBarVertStop = uint16(size(imageRepBright,1)*(1-numFractionScaleBarSpacing));
    numScaleBarVertStart = uint16(numScaleBarVertStop - numFractScaleBarWidth*size(imageRepBright,1));
    imageRepOut(numScaleBarVertStart:numScaleBarVertStop, numScaleBarHorStart:numScaleBarHorStop, :) = 255;
    
    %load the cropped representative image
    [ imRepCropRot, ~ ] = rotate_image( arrayZoomedRepImageRotate(iOutputProtein), double(imageRepRot), [1; 1] );
    imRepCropRot = uint8( imRepCropRot );
    imRepCropRot = imRepCropRot( arrayZoomedRepImageCrop(iOutputProtein,1):arrayZoomedRepImageCrop(iOutputProtein,2), ...
                                 arrayZoomedRepImageCrop(iOutputProtein,3):arrayZoomedRepImageCrop(iOutputProtein,4) );
    imRepCropRotBright = alterPixelIntensityForDisplay( imRepCropRot, numMinBrightnessThresh, numNonLinearTransCoefficient, numMaxBrightnessScaleCoefficient );
    imRepCropOut = cat( 3, imRepCropRotBright*0, imRepCropRotBright, imRepCropRotBright*0 );
    
    %draw a scale bar at the bottom right
    numScaleBarHorStop = uint16(size(imRepCropRotBright,2)*(1-numFractionScaleBarSpacing));
    numScaleBarHorStart = uint16(numScaleBarHorStop - arrayNumPixelsIn10um(iOutputProtein));
    numScaleBarVertStop = uint16(size(imRepCropRotBright,1)*(1-numFractionScaleBarSpacing));
    numScaleBarVertStart = uint16(numScaleBarVertStop - numFractScaleBarWidth*size(imRepCropRotBright,1));
    imRepCropOut(numScaleBarVertStart:numScaleBarVertStop, numScaleBarHorStart:numScaleBarHorStop, :) = 255;
    
    %load the 'other image', as specified at the start of this script
    structOtherImageOne = struct('dimensionality', {}, 'arrayXCoOrds', {}, 'arrayYCoOrds', {}, 'arrayPixelInt', {}, 'fv1', {}, 'fv2', {}, 'fv3', {});
    structOtherImageTwo = struct('dimensionality', {}, 'arrayXCoOrds', {}, 'arrayYCoOrds', {}, 'arrayPixelInt', {}, 'fv1', {}, 'fv2', {}, 'fv3', {});
    
    if (size(arrayOtherImageTypes{iOutputProtein},2) == 1),
    %display one image in the remaining space
        if arrayOtherImageTypes{iOutputProtein}(1) == 2,
            %plot the signal intensity as a surface
            structOtherImageOne(iOutputProtein).dimensionality = 2;
            
            %imageSurf = imrotate(imageRepRot, arrayOtherImageRot(iOutputProtein, 1));
            [imageSurf,~] = rotate_image( arrayOtherImageRot(iOutputProtein, 1), double(imageRepRot), [1; 1] );
            imageSurf = uint8(imageSurf);
            
            structOtherImageOne(iOutputProtein).arrayXCoOrds = [ arrayOtherImageCoOrds{iOutputProtein,1}(1):arrayOtherImageCoOrds{iOutputProtein,1}(2) ];
            structOtherImageOne(iOutputProtein).arrayYCoOrds = [ arrayOtherImageCoOrds{iOutputProtein,1}(3):arrayOtherImageCoOrds{iOutputProtein,1}(4) ];
            structOtherImageOne(iOutputProtein).arrayPixelInt = [ imageSurf( arrayOtherImageCoOrds{iOutputProtein,1}(3):arrayOtherImageCoOrds{iOutputProtein,1}(4), arrayOtherImageCoOrds{iOutputProtein,1}(1):arrayOtherImageCoOrds{iOutputProtein,1}(2) ) ];
            
        elseif arrayOtherImageTypes{iOutputProtein}(1) == 3,
            %plot the signal intensity as a volume
            structOtherImageOne(iOutputProtein).dimensionality = 3;
            numIsoSurfGroups = size(arrayIsoSurfGroups{iOutputProtein},1);
            
            %load the input image stack and isolate the target sub-volume
            arrayBSlashIndex = strfind(stringInputImagePath, '\');
            numLastBSlashIndex = arrayBSlashIndex(end);
            stringStackName = stringInputImagePath(numLastBSlashIndex+1:end);
            numZPosIndex = strfind(stringStackName, '_z');
            stringStackName(numZPosIndex:numZPosIndex+4) = '_z00*';
            numChIndex = strfind(stringStackName, '_ch');
            stringStackName(numChIndex:numChIndex+4) = '_ch0#';
            imageStackInput = loadImageStack( stringInputImagePath(1:numLastBSlashIndex), stringStackName, arrayImagesInStack(iOutputProtein), 0 );
            %rotate the image stack as required
            [imageTempInputRot, ~] = rotate_image( (arrayRepImageRotate(iOutputProtein)+arrayOtherImageRot(iOutputProtein, 1)), double(imageStackInput(:,:,1)), [1; 1] );
            imageStackInputRot = zeros(size(imageTempInputRot,1),size(imageTempInputRot,2),'uint8');
            for iImage = 1:arrayImagesInStack(iOutputProtein)
                [imageStackInputRot(:,:,iImage), ~] = rotate_image( (arrayRepImageRotate(iOutputProtein)+arrayOtherImageRot(iOutputProtein, 1)), double(imageStackInput(:,:,iImage)), [1; 1] );
            end
            imageStackSubVolume = imageStackInputRot(  arrayOtherImageCoOrds{iOutputProtein,1}(3):arrayOtherImageCoOrds{iOutputProtein,1}(4), ...
                                                       arrayOtherImageCoOrds{iOutputProtein,1}(1):arrayOtherImageCoOrds{iOutputProtein,1}(2), ...
                                                       arrayOtherImageCoOrds{iOutputProtein,1}(5):arrayOtherImageCoOrds{iOutputProtein,1}(6) );
            clear imageStackInput imageStackInputRot; %large array - clear early to free memory
            
            imageIndexedSubVolume = imageStackSubVolume*0;
            for iIsoSurfGroup = 1:numIsoSurfGroups,
                arrayVoxelsInGroup = find( imageStackSubVolume >= arrayIsoSurfGroups{iOutputProtein}{iIsoSurfGroup}(1) & ...
                                           imageStackSubVolume <= arrayIsoSurfGroups{iOutputProtein}{iIsoSurfGroup}(2) );
                imageIndexedSubVolume(arrayVoxelsInGroup) = iIsoSurfGroup;
                  
            end
            
            structOtherImageOne(iOutputProtein).fv1 = isosurface(imageIndexedSubVolume, 0);
            structOtherImageOne(iOutputProtein).fv2 = isosurface(imageIndexedSubVolume, 1);
            structOtherImageOne(iOutputProtein).fv3 = isosurface(imageIndexedSubVolume, 2);
            
        end
        
    elseif (size(arrayOtherImageTypes{iOutputProtein},2) == 2),
    %display two images in the remaining space
        %first 'other image'
        if arrayOtherImageTypes{iOutputProtein}(1) == 2,
            %plot the signal intensity as a surface
            structOtherImageOne(iOutputProtein).dimensionality = 2;
            
            %imageSurf = imrotate(imageRepRot, arrayOtherImageRot(iOutputProtein, 1));
            [imageSurf,JUNKarrayRefPoints] = rotate_image( arrayOtherImageRot(iOutputProtein, 1), double(imageRepRot), [1; 1] );
            imageSurf = uint8(imageSurf);
            
            structOtherImageOne(iOutputProtein).arrayXCoOrds = [ arrayOtherImageCoOrds{iOutputProtein,1}(1):arrayOtherImageCoOrds{iOutputProtein,1}(2) ];
            structOtherImageOne(iOutputProtein).arrayYCoOrds = [ arrayOtherImageCoOrds{iOutputProtein,1}(3):arrayOtherImageCoOrds{iOutputProtein,1}(4) ];
            structOtherImageOne(iOutputProtein).arrayPixelInt = [ imageSurf( arrayOtherImageCoOrds{iOutputProtein,1}(3):arrayOtherImageCoOrds{iOutputProtein,1}(4), arrayOtherImageCoOrds{iOutputProtein,1}(1):arrayOtherImageCoOrds{iOutputProtein,1}(2) ) ];
            
        elseif arrayOtherImageTypes{iOutputProtein}(1) == 3,
            %plot the signal intensity as a volume
            structOtherImageOne(iOutputProtein).dimensionality = 3;
            numIsoSurfGroups = size(arrayIsoSurfGroups{iOutputProtein},1);
            
            %load the input image stack and isolate the target sub-volume
            arrayBSlashIndex = strfind(stringInputImagePath, '\');
            numLastBSlashIndex = arrayBSlashIndex(end);
            stringStackName = stringInputImagePath(numLastBSlashIndex+1:end);
            numZPosIndex = strfind(stringStackName, '_z');
            stringStackName(numZPosIndex:numZPosIndex+4) = '_z00*';
            numChIndex = strfind(stringStackName, '_ch');
            stringStackName(numChIndex:numChIndex+4) = '_ch0#';
            imageStackInput = loadImageStack( stringInputImagePath(1:numLastBSlashIndex), stringStackName, arrayImagesInStack(iOutputProtein), 0 );
            %rotate the image stack as required
            [imageTempInputRot, JUNKarrayRefPoints] = rotate_image( (arrayRepImageRotate(iOutputProtein)+arrayOtherImageRot(iOutputProtein, 1)), double(imageStackInput(:,:,1)), [1; 1] );
            imageStackInputRot = zeros(size(imageTempInputRot,1),size(imageTempInputRot,2),'uint8');
            for iImage = 1:arrayImagesInStack(iOutputProtein)
                [imageStackInputRot(:,:,iImage), JUNKarrayRefPoints] = rotate_image( (arrayRepImageRotate(iOutputProtein)+arrayOtherImageRot(iOutputProtein, 1)), double(imageStackInput(:,:,iImage)), [1; 1] );
            end
            imageStackSubVolume = imageStackInputRot(  arrayOtherImageCoOrds{iOutputProtein,1}(3):arrayOtherImageCoOrds{iOutputProtein,1}(4), ...
                                                       arrayOtherImageCoOrds{iOutputProtein,1}(1):arrayOtherImageCoOrds{iOutputProtein,1}(2), ...
                                                       arrayOtherImageCoOrds{iOutputProtein,1}(5):arrayOtherImageCoOrds{iOutputProtein,1}(6) );
            clear imageStackInput imageStackInputRot; %large array - clear early to free memory
            
            imageIndexedSubVolume = imageStackSubVolume*0;
            for iIsoSurfGroup = 1:numIsoSurfGroups,
                arrayVoxelsInGroup = find( imageStackSubVolume >= arrayIsoSurfGroups{iOutputProtein}{iIsoSurfGroup}(1) & ...
                                           imageStackSubVolume <= arrayIsoSurfGroups{iOutputProtein}{iIsoSurfGroup}(2) );
                imageIndexedSubVolume(arrayVoxelsInGroup) = iIsoSurfGroup;
                  
            end
            
            structOtherImageOne(iOutputProtein).fv1 = isosurface(imageIndexedSubVolume, 0);
            structOtherImageOne(iOutputProtein).fv2 = isosurface(imageIndexedSubVolume, 1);
            structOtherImageOne(iOutputProtein).fv3 = isosurface(imageIndexedSubVolume, 2);
            
        end
        %second 'other image'
        if arrayOtherImageTypes{iOutputProtein}(2) == 2,
            %plot the signal intensity as a surface
            structOtherImageTwo(iOutputProtein).dimensionality = 2;
            
            %imageSurf = imrotate(imageRepRot, arrayOtherImageRot(iOutputProtein, 1));
            [imageSurf,JUNKarrayRefPoints] = rotate_image( arrayOtherImageRot(iOutputProtein, 2), double(imageRepRot), [1; 1] );
            imageSurf = uint8(imageSurf);
            
            structOtherImageTwo(iOutputProtein).arrayXCoOrds = [ arrayOtherImageCoOrds{iOutputProtein,2}(1):arrayOtherImageCoOrds{iOutputProtein,2}(2) ];
            structOtherImageTwo(iOutputProtein).arrayYCoOrds = [ arrayOtherImageCoOrds{iOutputProtein,2}(3):arrayOtherImageCoOrds{iOutputProtein,2}(4) ];
            structOtherImageTwo(iOutputProtein).arrayPixelInt = [ imageSurf( arrayOtherImageCoOrds{iOutputProtein,2}(3):arrayOtherImageCoOrds{iOutputProtein,2}(4), arrayOtherImageCoOrds{iOutputProtein,2}(1):arrayOtherImageCoOrds{iOutputProtein,2}(2) ) ];
            
        elseif arrayOtherImageTypes{iOutputProtein}(2) == 3,
            %plot the signal intensity as a volume
            structOtherImageTwo(iOutputProtein).dimensionality = 3;
            numIsoSurfGroups = size(arrayIsoSurfGroups{iOutputProtein},1);
            
            %load the input image stack and isolate the target sub-volume
            arrayBSlashIndex = strfind(stringInputImagePath, '\');
            numLastBSlashIndex = arrayBSlashIndex(end);
            stringStackName = stringInputImagePath(numLastBSlashIndex+1:end);
            numZPosIndex = strfind(stringStackName, '_z');
            stringStackName(numZPosIndex:numZPosIndex+4) = '_z00*';
            numChIndex = strfind(stringStackName, '_ch');
            stringStackName(numChIndex:numChIndex+4) = '_ch0#';
            imageStackInput = loadImageStack( stringInputImagePath(1:numLastBSlashIndex), stringStackName, arrayImagesInStack(iOutputProtein), 0 );
            %rotate the image stack as required
            [imageTempInputRot, JUNKarrayRefPoints] = rotate_image( (arrayRepImageRotate(iOutputProtein)+arrayOtherImageRot(iOutputProtein, 2)), double(imageStackInput(:,:,1)), [1; 1] );
            imageStackInputRot = zeros(size(imageTempInputRot,1),size(imageTempInputRot,2),'uint8');
            for iImage = 1:arrayImagesInStack(iOutputProtein)
                [imageStackInputRot(:,:,iImage), JUNKarrayRefPoints] = rotate_image( (arrayRepImageRotate(iOutputProtein)+arrayOtherImageRot(iOutputProtein, 2)), double(imageStackInput(:,:,iImage)), [1; 1] );
            end
            imageStackSubVolume = imageStackInputRot(  arrayOtherImageCoOrds{iOutputProtein,2}(3):arrayOtherImageCoOrds{iOutputProtein,2}(4), ...
                                                       arrayOtherImageCoOrds{iOutputProtein,2}(1):arrayOtherImageCoOrds{iOutputProtein,2}(2), ...
                                                       arrayOtherImageCoOrds{iOutputProtein,2}(5):arrayOtherImageCoOrds{iOutputProtein,2}(6) );
            clear imageStackInput imageStackInputRot; %large array - clear early to free memory
            
            imageIndexedSubVolume = imageStackSubVolume*0;
            for iIsoSurfGroup = 1:numIsoSurfGroups,
                arrayVoxelsInGroup = find( imageStackSubVolume >= arrayIsoSurfGroups{iOutputProtein}{iIsoSurfGroup}(1) & ...
                                           imageStackSubVolume <= arrayIsoSurfGroups{iOutputProtein}{iIsoSurfGroup}(2) );
                imageIndexedSubVolume(arrayVoxelsInGroup) = iIsoSurfGroup;
                  
            end
            
            structOtherImageTwo(iOutputProtein).fv1 = isosurface(imageIndexedSubVolume, 0);
            structOtherImageTwo(iOutputProtein).fv2 = isosurface(imageIndexedSubVolume, 1);
            structOtherImageTwo(iOutputProtein).fv3 = isosurface(imageIndexedSubVolume, 2);
            
        end
    end
    
    %calculate the position of some reference points to determine the
    %transform from the representative image to the zoomed image
    arrayOriginalImagePoints = [ 0001, 0001;  %x y      % origin
                                 0001, 1024 ];          % bottom left
    if arrayZoomedRepImageRotate(iOutputProtein) < 0,
        arrayRotatedImagePoints = [ 1+numImageSizePixels*sind(180+arrayZoomedRepImageRotate(iOutputProtein)), 0001;         % origin
                                    0001, 1+numImageSizePixels*cosd(arrayZoomedRepImageRotate(iOutputProtein))  ];      % bottom left
        %convert to integers and back to double
        arrayRotatedImagePoints = double(uint16(arrayRotatedImagePoints));
    elseif arrayZoomedRepImageRotate(iOutputProtein) == 0,
        %there is no transform
        arrayRotatedImagePoints = [ 0001, 0001;  %x y      % origin
                                     0001, 1024 ];          % bottom left
    elseif arrayZoomedRepImageRotate(iOutputProtein) > 0,
        arrayRotatedImagePoints = [ 0001, 1+numImageSizePixels*sind(arrayZoomedRepImageRotate(iOutputProtein));         % origin
                                    1+numImageSizePixels*sind(arrayZoomedRepImageRotate(iOutputProtein)), 1+numImageSizePixels*(sind(arrayZoomedRepImageRotate(iOutputProtein))+cosd(arrayZoomedRepImageRotate(iOutputProtein)))  ];      % bottom left
        %convert to integers and back to double
        arrayRotatedImagePoints = double(uint16(arrayRotatedImagePoints));
    else
        disp('the rotation angle is not right');
    end
    structRotationTransform = cp2tform(arrayOriginalImagePoints, arrayRotatedImagePoints, 'nonreflective similarity');
    
    %calculate the zoomed region co-ordinates on the original image
    arrayZoomedImageXCoOrd = [ arrayZoomedRepImageCrop(iOutputProtein,3);
                               arrayZoomedRepImageCrop(iOutputProtein,4);
                               arrayZoomedRepImageCrop(iOutputProtein,4);
                               arrayZoomedRepImageCrop(iOutputProtein,3)    ];
    arrayZoomedImageYCoOrd = [ arrayZoomedRepImageCrop(iOutputProtein,1);
                               arrayZoomedRepImageCrop(iOutputProtein,1);
                               arrayZoomedRepImageCrop(iOutputProtein,2);
                               arrayZoomedRepImageCrop(iOutputProtein,2)    ];
    [  arrayZoomedImgOrigXCoOrd, arrayZoomedImgOrigYCoOrd  ] = ...
        tforminv(structRotationTransform, arrayZoomedImageXCoOrd, arrayZoomedImageYCoOrd);
                
    %calculate a similar back-transform of the ROI for the 'other images'
    arrayOriginalImagePoints = [ 0001, 0001;  %x y      % origin
                                 0001, 1024 ];          % bottom left
    if arrayOtherImageRot(iOutputProtein, 1) < 0,
        arrayRotatedImagePoints = [ 1+numImageSizePixels*sind(180+arrayOtherImageRot(iOutputProtein, 1)), 0001;         % origin
                                    0001, 1+numImageSizePixels*cosd(arrayOtherImageRot(iOutputProtein, 1))  ];      % bottom left
        %convert to integers and back to double
        arrayRotatedImagePoints = double(uint16(arrayRotatedImagePoints));
    elseif arrayOtherImageRot(iOutputProtein, 1) == 0,
        %there is no transform
        arrayRotatedImagePoints = [ 0001, 0001;  %x y      % origin
                                     0001, 1024 ];          % bottom left
    elseif arrayOtherImageRot(iOutputProtein, 1) > 0,
        arrayRotatedImagePoints = [ 0001, 1+numImageSizePixels*sind(arrayOtherImageRot(iOutputProtein, 1));         % origin
                                    1+numImageSizePixels*sind(arrayOtherImageRot(iOutputProtein, 1)), 1+numImageSizePixels*(sind(arrayOtherImageRot(iOutputProtein, 1))+cosd(arrayOtherImageRot(iOutputProtein, 1)))  ];      % bottom left
        %convert to integers and back to double
        arrayRotatedImagePoints = double(uint16(arrayRotatedImagePoints));
    else
        disp('the rotation angle is not right');
    end
    structRotationTransform = cp2tform(arrayOriginalImagePoints, arrayRotatedImagePoints, 'nonreflective similarity');

    %calculate the zoomed region co-ordinates on the original image
    arrayOtherImageXCoOrd = [ arrayOtherImageCoOrds{iOutputProtein,1}(1);
                               arrayOtherImageCoOrds{iOutputProtein,1}(2);
                               arrayOtherImageCoOrds{iOutputProtein,1}(2);
                               arrayOtherImageCoOrds{iOutputProtein,1}(1)    ];
    arrayOtherImageYCoOrd = [ arrayOtherImageCoOrds{iOutputProtein,1}(3);
                               arrayOtherImageCoOrds{iOutputProtein,1}(3);
                               arrayOtherImageCoOrds{iOutputProtein,1}(4);
                               arrayOtherImageCoOrds{iOutputProtein,1}(4)    ];
    [  arrayOtherImgOneOrigXCoOrd, arrayOtherImgOneOrigYCoOrd  ] = ...
        tforminv(structRotationTransform, arrayOtherImageXCoOrd, arrayOtherImageYCoOrd);
    
    %calculate the approximate d_norm values for the zoomed region
    if size(arrayRepImageZ{iOutputProtein},2) > 1,
        numZSlicesInRepZ = size(arrayRepImageZ{iOutputProtein},2);
        numZToSearch = arrayRepImageZ{iOutputProtein}(ceil(numZSlicesInRepZ/2));
    else
        numZToSearch = arrayRepImageZ{iOutputProtein};
    end
    if isempty(  find( arrayZSlicesSampled == numZToSearch, 1 )  ),
        %find the closest sampled z-slice, to the one displayed
        [ JUNKnumMin, arrayMinIndex ] = min(abs(arrayZSlicesSampled - numZToSearch));
        numZSampleForDNormEst = arrayZSlicesSampled(arrayMinIndex);
    else
        numZSampleForDNormEst = numZToSearch;
    end
    
    [ CellImageStack, CellImageLayers, CellImageBasalLamina, CellImageBoundaryOne, CellImageBoundaryTwo, CellImageOuterBoundary ] = ...
        loadImageStackAsSparse3DCellArrays( numZSampleForDNormEst, stringSpecificImageDataFolder, stringRepImageDataFolder);
    arrayBorderLocations = struct('SmpCent', {}, 'ZPosition', {});
    numXCoOrd1 = arrayOtherImageXCoOrd(4);
    numXCoOrd2 = arrayOtherImageXCoOrd(3);
    numYCoOrd1 = arrayOtherImageYCoOrd(4);
    numYCoOrd2 = arrayOtherImageYCoOrd(3);
    [ JUNKimageRepRotInv, arrayRefPoints ] = rotate_image( -arrayRepImageRotate(iOutputProtein), double(imageRepRot), [numXCoOrd1, numXCoOrd2; numYCoOrd1, numYCoOrd2] );
    arrayBorderLocations(1).SmpCent = [ arrayRefPoints(1,1) arrayRefPoints(2,1)  numZSampleForDNormEst ];
    arrayBorderLocations(2).SmpCent = [ arrayRefPoints(1,2) arrayRefPoints(2,2)  numZSampleForDNormEst ];
    arrayBorderLocations(1).ZPosition = numZSampleForDNormEst;
    arrayBorderLocations(2).ZPosition = numZSampleForDNormEst;
    arraySamplingKernel = [ 1 ];
    
    arrayOutSampleAnalysis = loadSampAnalysis(stringRepImageCytoDataFolder, '.mat');
    
    %size of the output figure
    numFigHeight = (55+(50*numLocalisationsForProtein));
    if (iOutputProtein >= 7) && (iOutputProtein <= 10),
        %MEK, pMEK, ERK or pERK; plot the cyto:nuc signal ratio in the
        % final panel
        numFigHeight = (55+(50*(numLocalisationsForProtein+1)));
    end
    numFigWidth = 150;
    arrayFigPosition = [ 200, 300, numFigScaleMult*numFigWidth, numFigScaleMult*numFigHeight ];
    
    
    %representative image position
    arrayRepImagePosition = [ 5/numFigWidth,  (numFigHeight-48)/numFigHeight,  45/numFigWidth, 45/numFigHeight ];
    arrayZoomedRepImagePosition = [ 60/numFigWidth,  (numFigHeight-48)/numFigHeight,  80/numFigWidth, 45/numFigHeight ];

    
    %overlaid lowess curve positions
    arrayPlot1Position = [ 80/numFigWidth,  (numFigHeight-93)/numFigHeight,  60/numFigWidth,  35/numFigHeight ];
    arrayPlot2Position = [ 80/numFigWidth,  (numFigHeight-143)/numFigHeight,  60/numFigWidth,  35/numFigHeight ];
    arrayPlot3Position = [ 80/numFigWidth,  (numFigHeight-193)/numFigHeight,  60/numFigWidth,  35/numFigHeight ];
    
    
    %'other plot' positions (volume and surface plots)
    arrayOneLocOtherImagePosition = [ 5/numFigWidth, (numFigHeight-98)/numFigHeight, 45/numFigWidth, 45/numFigHeight ];
    
    arrayTwoLocOtherImagePositionSingle = [ 5/numFigWidth, (numFigHeight-103)/numFigHeight, 50/numFigWidth, 50/numFigHeight ];
    arrayTwoLocOtherImagePositionDoubleOne = [ 5/numFigWidth, (numFigHeight-98)/numFigHeight, 45/numFigWidth, 45/numFigHeight ];
    arrayTwoLocOtherImagePositionDoubleTwo = [ 5/numFigWidth, (numFigHeight-148)/numFigHeight, 45/numFigWidth, 45/numFigHeight ];
    
    arrayThreeLocOtherImagePositionOne = [ 5/numFigWidth, (numFigHeight-103)/numFigHeight, 50/numFigWidth, 50/numFigHeight ];
    arrayThreeLocOtherImagePositionTwo = [ 5/numFigWidth, (numFigHeight-168)/numFigHeight, 50/numFigWidth, 50/numFigHeight ];
    
    %initialise the output figure
    figOut = figure;
    set(figOut, 'Position', arrayFigPosition);
    set(figOut, 'PaperUnits', 'centimeters', 'PaperSize', [ 15, 5.5+(4.5*numLocalisationsForProtein) ], 'PaperPosition', [0, 0, 15, 5.5+(4.5*numLocalisationsForProtein)] );
    if (iOutputProtein >= 7) && (iOutputProtein <= 10),
        %MEK, pMEK, ERK or pERK; plot the cyto:nuc signal ratio in the
        % final panel
        set(figOut, 'PaperUnits', 'centimeters', 'PaperSize', [ 15, 5.5+(4.5*(numLocalisationsForProtein+1)) ], 'PaperPosition', [0, 0, 15, 5.5+(4.5*(numLocalisationsForProtein+1))] );
    end    
    
    %display the unzoomed representative image
    handleRepImage = subplot('Position', arrayRepImagePosition);
    image(imageRepOut);
    set(handleRepImage, 'XTick', [], 'YTick', []);
    text(size(imageRepOut,1)/23,size(imageRepOut,2)/11,'A','FontSize',numSubFigLabelFontSize,'Color','w', 'FontWeight', 'bold');
    %overlay the borders of the zoomed region, on the unzoomed image
    line([arrayZoomedImgOrigXCoOrd(1), arrayZoomedImgOrigXCoOrd(2)], [arrayZoomedImgOrigYCoOrd(1), arrayZoomedImgOrigYCoOrd(2)], 'LineWidth', numOverlayBorderWidth, 'Color', 'w', 'LineStyle', ':');
    line([arrayZoomedImgOrigXCoOrd(2), arrayZoomedImgOrigXCoOrd(3)], [arrayZoomedImgOrigYCoOrd(2), arrayZoomedImgOrigYCoOrd(3)], 'LineWidth', numOverlayBorderWidth, 'Color', 'w', 'LineStyle', ':');
    line([arrayZoomedImgOrigXCoOrd(3), arrayZoomedImgOrigXCoOrd(4)], [arrayZoomedImgOrigYCoOrd(3), arrayZoomedImgOrigYCoOrd(4)], 'LineWidth', numOverlayBorderWidth, 'Color', 'w', 'LineStyle', ':');
    line([arrayZoomedImgOrigXCoOrd(4), arrayZoomedImgOrigXCoOrd(1)], [arrayZoomedImgOrigYCoOrd(4), arrayZoomedImgOrigYCoOrd(1)], 'LineWidth', numOverlayBorderWidth, 'Color', 'w', 'LineStyle', ':');
    %overlay the borders of the 'other image' region, on the unzoomed image
    line([arrayOtherImgOneOrigXCoOrd(1), arrayOtherImgOneOrigXCoOrd(2)], [arrayOtherImgOneOrigYCoOrd(1), arrayOtherImgOneOrigYCoOrd(2)], 'LineWidth', numOverlayBorderWidth, 'Color', [1 0.5 0], 'LineStyle', ':');
    line([arrayOtherImgOneOrigXCoOrd(2), arrayOtherImgOneOrigXCoOrd(3)], [arrayOtherImgOneOrigYCoOrd(2), arrayOtherImgOneOrigYCoOrd(3)], 'LineWidth', numOverlayBorderWidth, 'Color', [1 0.5 0], 'LineStyle', ':');
    line([arrayOtherImgOneOrigXCoOrd(3), arrayOtherImgOneOrigXCoOrd(4)], [arrayOtherImgOneOrigYCoOrd(3), arrayOtherImgOneOrigYCoOrd(4)], 'LineWidth', numOverlayBorderWidth, 'Color', [1 0.5 0], 'LineStyle', '-');
    line([arrayOtherImgOneOrigXCoOrd(4), arrayOtherImgOneOrigXCoOrd(1)], [arrayOtherImgOneOrigYCoOrd(4), arrayOtherImgOneOrigYCoOrd(1)], 'LineWidth', numOverlayBorderWidth, 'Color', [1 0.5 0], 'LineStyle', ':');
    
    
    
    %display the zoomed representative image
    handleRepZoomedImage = subplot('Position', arrayZoomedRepImagePosition);
    image(imRepCropOut);
    set(handleRepZoomedImage, 'XTick', [], 'YTick', []);
    text(size(imRepCropOut,1)/30,size(imRepCropOut,2)/19,'B','FontSize',numSubFigLabelFontSize,'Color','w', 'FontWeight', 'bold');
    %insert the protein name
    numTextXPos = size(imRepCropOut,2)*0.9;
    numTextYPos = size(imRepCropOut,1)*0.13;
    text(numTextXPos, numTextYPos, arrayProteinLateXString{iOutputProtein}, 'Interpreter', 'latex', 'FontWeight', 'bold', 'HorizontalAlignment', 'right','FontSize',(numSubFigLabelFontSize*0.85),'Color','w');
    
 % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% Plot the 2D surface or 3D volume
 % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    %display the 'other image'
    if numLocalisationsForProtein == 1,
        %plot a single 'other image'
        handleOtherImage = subplot('Position', arrayOneLocOtherImagePosition);
        
        if structOtherImageOne(iOutputProtein).dimensionality == 2,
            %surface plot
            arrayXRange = structOtherImageOne(iOutputProtein).arrayXCoOrds;
            arrayYRange = structOtherImageOne(iOutputProtein).arrayYCoOrds;
            arrayXRangeInterp = arrayXRange(1):(1/double(arrayNumSurfaceInterp(iOutputProtein,1))):arrayXRange(end);
            arrayYRangeInterp = arrayYRange(1):(1/double(arrayNumSurfaceInterp(iOutputProtein,1))):arrayYRange(end);
            imageSubArea = structOtherImageOne(iOutputProtein).arrayPixelInt;
            imageSubAreaInterp = interp2(arrayXRange, arrayYRange, ...
                                         double(imageSubArea(:,:)), ...
                                         arrayXRangeInterp, arrayYRangeInterp', ...
                                         'cubic');
            
            arraySubAreaInterpBelowZeroObsIndices = find(imageSubAreaInterp < 0.01);
            imageSubAreaInterp(arraySubAreaInterpBelowZeroObsIndices) = 0.01;

            if isdeployed,
                %don't use alpha data or lighting, because this can crash 
                % the virtual box (due to OpenGL rendering issues)
                surf(arrayXRangeInterp, arrayYRangeInterp, imageSubAreaInterp, ...
                    'FaceColor','interp',  'EdgeColor','none');
            else
                %plot the surface using alpha transparency
                arrayAlphaData = log(imageSubAreaInterp);
                surf(arrayXRangeInterp, arrayYRangeInterp, imageSubAreaInterp, ...
                    'FaceColor','interp',  'EdgeColor','none',  'FaceAlpha','interp',...
                    'AlphaDataMapping','scaled',  'AlphaData',arrayAlphaData, ...
                    'FaceLighting','phong', 'AmbientStrength',0.8);
                camlight;
            end
            axis tight;
            axis ij;
            axis equal;
            colormap(jet);
            set(gca, 'DataAspectRatio',[1 1 4/double(arrayNumSurfaceInterp(iOutputProtein,1))]);
            view(arrayOtherImageView{iOutputProtein,1}(1), arrayOtherImageView{iOutputProtein,1}(2));
            set(gca, 'YTick', [], 'YTickLabel', [], 'ZTick', [], 'ZTickLabel', []);
            set(gca, 'XTick', [arrayXRangeInterp(1) arrayXRangeInterp(end)], 'XTickLabel', [{['~' num2str(arrayOutSampleAnalysis(1).NormDist, '%3.2f')]}; {['~' num2str(arrayOutSampleAnalysis(2).NormDist, '%3.2f')]}]);
            text((arrayXRange(1) + arrayXRange(end))/2, arrayXRange(end)*1.35, 0, {'Approx.';'d_n_o_r_m'}, 'HorizontalAlignment', 'center');
            zlabel('Pix. Int.');
            line([min(arrayXRangeInterp) max(arrayXRangeInterp)], [max(arrayYRangeInterp) max(arrayYRangeInterp)], [0 0], 'LineWidth', numOverlayBorderWidth, 'Color', [1 0.5 0], 'LineStyle', '-');
            text(min(arrayXRangeInterp),min(arrayYRangeInterp),max(imageSubAreaInterp(:))*1.5,'C','FontSize',numSubFigLabelFontSize,'Color','k', 'FontWeight', 'bold');
            
        elseif structOtherImageOne(iOutputProtein).dimensionality == 3,
            %volume plot
            hold on;
            if isdeployed,
                %don't use alpha data or lighting, because this can crash 
                % the virtual box (due to OpenGL rendering issues)
                plotPatch1 = patch(structOtherImageOne(iOutputProtein).fv1, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none' );

                plotPatch2 = patch(structOtherImageOne(iOutputProtein).fv2, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none' );

                plotPatch3 = patch(structOtherImageOne(iOutputProtein).fv3, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none' );
            else
                %plot the surface using alpha transparency
                plotPatch1 = patch(structOtherImageOne(iOutputProtein).fv1, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none', ...
                    'FaceAlpha', 0.05, ...
                    'FaceLighting','phong', ...
                    'AmbientStrength',0.1);

                plotPatch2 = patch(structOtherImageOne(iOutputProtein).fv2, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none', ...
                    'FaceAlpha', 0.5, ...
                    'FaceLighting','phong', ...
                    'AmbientStrength',0.5);

                plotPatch3 = patch(structOtherImageOne(iOutputProtein).fv3, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none', ...
                    'FaceAlpha', 0.7, ...
                    'FaceLighting','phong', ...
                    'AmbientStrength', 0.7);
                
                %turn on the cam light
                camlight;
            end
            
            line( [0 arrayOtherImageCoOrds{iOutputProtein,1}(2)-arrayOtherImageCoOrds{iOutputProtein,1}(1)], ...
                  [arrayOtherImageCoOrds{iOutputProtein,1}(4)-arrayOtherImageCoOrds{iOutputProtein,1}(3) arrayOtherImageCoOrds{iOutputProtein,1}(4)-arrayOtherImageCoOrds{iOutputProtein,1}(3)], ...
                  [0 0], 'LineWidth', numOverlayBorderWidth, 'Color', [1 0.5 0], 'LineStyle', '-');
            
            box on;
            axis tight;
            axis ij;
            view(arrayOtherImageView{3,1});
            axis equal;
            set(gca, 'Projection', 'Perspective');
            set(gca, 'YTick', [], 'YTickLabel', [], 'ZTick', [], 'ZTickLabel', []);
            set(gca, 'XTick', [arrayXRangeInterp(1) arrayXRangeInterp(end)], 'XTickLabel', [{['~' num2str(arrayOutSampleAnalysis(1).NormDist, '%3.2f')]}; {['~' num2str(arrayOutSampleAnalysis(2).NormDist, '%3.2f')]}]);
            xlabel('Approx. d_n_o_r_m');
            set(gca, 'DataAspectRatio', [1 1 0.6]);
            hold off;
            text(0,0,(arrayOtherImageCoOrds{iOutputProtein,1}(6)-arrayOtherImageCoOrds{iOutputProtein,1}(5))*1.5,'C','FontSize',numSubFigLabelFontSize,'Color','k', 'FontWeight', 'bold');
            
        end
        
    elseif (numLocalisationsForProtein == 2) && (length(arrayOtherImageTypes{iOutputProtein}) == 1),
        %plot a single 'other image'
        handleOtherImage = subplot('Position', arrayTwoLocOtherImagePositionSingle);
        if structOtherImageOne(iOutputProtein).dimensionality == 2,
            %surface plot
            arrayXRange = structOtherImageOne(iOutputProtein).arrayXCoOrds;
            arrayYRange = structOtherImageOne(iOutputProtein).arrayYCoOrds;
            arrayXRangeInterp = arrayXRange(1):(1/double(arrayNumSurfaceInterp(iOutputProtein,1))):arrayXRange(end);
            arrayYRangeInterp = arrayYRange(1):(1/double(arrayNumSurfaceInterp(iOutputProtein,1))):arrayYRange(end);
            imageSubArea = structOtherImageOne(iOutputProtein).arrayPixelInt;
            imageSubAreaInterp = interp2(arrayXRange, arrayYRange, ...
                                         double(imageSubArea(:,:)), ...
                                         arrayXRangeInterp, arrayYRangeInterp', ...
                                         'cubic');
            arraySubAreaInterpBelowZeroObsIndices = find(imageSubAreaInterp < 0.01);
            imageSubAreaInterp(arraySubAreaInterpBelowZeroObsIndices) = 0.01;
            if isdeployed,
                %don't use alpha data or lighting, because this can crash 
                % the virtual box (due to OpenGL rendering issues)
                surf(arrayXRangeInterp, arrayYRangeInterp, imageSubAreaInterp, ...
                    'FaceColor','interp',  'EdgeColor','none');
            else
                %plot the surface using alpha transparency
                arrayAlphaData = log(imageSubAreaInterp);
                surf(arrayXRangeInterp, arrayYRangeInterp, imageSubAreaInterp, ...
                    'FaceColor','interp',  'EdgeColor','none',  'FaceAlpha','interp',...
                    'AlphaDataMapping','scaled',  'AlphaData',arrayAlphaData, ...
                    'FaceLighting','phong', 'AmbientStrength',0.8);
                camlight;
            end
            axis tight;
            axis ij;
            axis equal;
            colormap(jet);
            set(gca, 'DataAspectRatio',[1 1 4/double(arrayNumSurfaceInterp(iOutputProtein,1))]);
            view(arrayOtherImageView{iOutputProtein,1}(1), arrayOtherImageView{iOutputProtein,1}(2));
            set(gca, 'YTick', [], 'YTickLabel', [], 'ZTick', [], 'ZTickLabel', []);
            set(gca, 'XTick', [arrayXRangeInterp(1) arrayXRangeInterp(end)], 'XTickLabel', [{['~' num2str(arrayOutSampleAnalysis(1).NormDist, '%3.2f')]}; {['~' num2str(arrayOutSampleAnalysis(2).NormDist, '%3.2f')]}]);
            zlabel('Pix. Int.');
            xlabel('Approx. d_n_o_r_m');
            line([min(arrayXRangeInterp) max(arrayXRangeInterp)], [max(arrayYRangeInterp) max(arrayYRangeInterp)], [0 0], 'LineWidth', numOverlayBorderWidth, 'Color', [1 0.5 0], 'LineStyle', '-');
            text(min(arrayXRangeInterp),min(arrayYRangeInterp),max(imageSubAreaInterp(:))*1.5,'C','FontSize',numSubFigLabelFontSize,'Color','k', 'FontWeight', 'bold');
            
        elseif structOtherImageOne(iOutputProtein).dimensionality == 3,
            %volume plot
            hold on;
            if isdeployed,
                %don't use alpha data or lighting, because this can crash 
                % the virtual box (due to OpenGL rendering issues)
                plotPatch1 = patch(structOtherImageOne(iOutputProtein).fv1, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none' );

                plotPatch2 = patch(structOtherImageOne(iOutputProtein).fv2, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none' );

                plotPatch3 = patch(structOtherImageOne(iOutputProtein).fv3, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none' );
            else
                %plot using alpha transparency
                plotPatch1 = patch(structOtherImageOne(iOutputProtein).fv1, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none', ...
                    'FaceAlpha', 0.05, ...
                    'FaceLighting','phong', ...
                    'AmbientStrength',0.1);

                plotPatch2 = patch(structOtherImageOne(iOutputProtein).fv2, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none', ...
                    'FaceAlpha', 0.5, ...
                    'FaceLighting','phong', ...
                    'AmbientStrength',0.5);

                plotPatch3 = patch(structOtherImageOne(iOutputProtein).fv3, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none', ...
                    'FaceAlpha', 0.7, ...
                    'FaceLighting','phong', ...
                    'AmbientStrength', 0.7);
                
                camlight;
            end
            
            line( [0 arrayOtherImageCoOrds{iOutputProtein,1}(2)-arrayOtherImageCoOrds{iOutputProtein,1}(1)], ...
                  [arrayOtherImageCoOrds{iOutputProtein,1}(4)-arrayOtherImageCoOrds{iOutputProtein,1}(3) arrayOtherImageCoOrds{iOutputProtein,1}(4)-arrayOtherImageCoOrds{iOutputProtein,1}(3)], ...
                  [0 0], 'LineWidth', numOverlayBorderWidth, 'Color', [1 0.5 0], 'LineStyle', '-');
            
            box on;
            axis tight;
            axis ij;
            view(arrayOtherImageView{3,1});
            axis equal;
            set(gca, 'Projection', 'Perspective');
            set(gca, 'YTick', [], 'YTickLabel', [], 'ZTick', [], 'ZTickLabel', []);
            set(gca, 'XTick', [0 arrayOtherImageCoOrds{iOutputProtein,1}(2)-arrayOtherImageCoOrds{iOutputProtein,1}(1)], 'XTickLabel', [{['~' num2str(arrayOutSampleAnalysis(1).NormDist, '%3.2f')]}; {['~' num2str(arrayOutSampleAnalysis(2).NormDist, '%3.2f')]}]);
            xlabel('Approx. d_n_o_r_m');
            set(gca, 'DataAspectRatio', [1 1 0.6]);
            hold off;
            text(0,0,(arrayOtherImageCoOrds{iOutputProtein,1}(6)-arrayOtherImageCoOrds{iOutputProtein,1}(5))*1.5,'C','FontSize',numSubFigLabelFontSize,'Color','k', 'FontWeight', 'bold');
            
        end
        
    elseif (numLocalisationsForProtein == 2) && (length(arrayOtherImageTypes{iOutputProtein}) == 2),
        %plot the first 'other image'
        handleOtherImage = subplot('Position', arrayTwoLocOtherImagePositionDoubleOne);
        if structOtherImageOne(iOutputProtein).dimensionality == 2,
            %surface plot
            arrayXRange = structOtherImageOne(iOutputProtein).arrayXCoOrds;
            arrayYRange = structOtherImageOne(iOutputProtein).arrayYCoOrds;
            arrayXRangeInterp = arrayXRange(1):(1/double(arrayNumSurfaceInterp(iOutputProtein,1))):arrayXRange(end);
            arrayYRangeInterp = arrayYRange(1):(1/double(arrayNumSurfaceInterp(iOutputProtein,1))):arrayYRange(end);
            imageSubArea = structOtherImageOne(iOutputProtein).arrayPixelInt;
            imageSubAreaInterp = interp2(arrayXRange, arrayYRange, ...
                                         double(imageSubArea(:,:)), ...
                                         arrayXRangeInterp, arrayYRangeInterp', ...
                                         'cubic');
            arraySubAreaInterpBelowZeroObsIndices = find(imageSubAreaInterp < 0.01);
            imageSubAreaInterp(arraySubAreaInterpBelowZeroObsIndices) = 0.01;
            
            if isdeployed,
                %don't use alpha data or lighting, because this can crash 
                % the virtual box (due to OpenGL rendering issues)
                surf(arrayXRangeInterp, arrayYRangeInterp, imageSubAreaInterp, ...
                    'FaceColor','interp',  'EdgeColor','none' );
            else
                
                arrayAlphaData = log(imageSubAreaInterp);

                surf(arrayXRangeInterp, arrayYRangeInterp, imageSubAreaInterp, ...
                    'FaceColor','interp',  'EdgeColor','none',  'FaceAlpha','interp',...
                    'AlphaDataMapping','scaled',  'AlphaData',arrayAlphaData, ...
                    'FaceLighting','phong', 'AmbientStrength',0.8);
                
                camlight;
            end
            axis tight;
            axis ij;
            axis equal;
            colormap(jet);
            set(gca, 'DataAspectRatio',[1 1 4/double(arrayNumSurfaceInterp(iOutputProtein,1))]);
            view(arrayOtherImageView{iOutputProtein,1}(1), arrayOtherImageView{iOutputProtein,1}(2));
            set(gca, 'YTick', [], 'YTickLabel', [], 'ZTick', [], 'ZTickLabel', []);
            set(gca, 'XTick', [arrayXRangeInterp(1) arrayXRangeInterp(end)], 'XTickLabel', [{['~' num2str(arrayOutSampleAnalysis(1).NormDist, '%3.2f')]}; {['~' num2str(arrayOutSampleAnalysis(2).NormDist, '%3.2f')]}]);
            zlabel('Pix. Int.');
            text((arrayXRange(1) + arrayXRange(end))/2, arrayYRange(end)*1.08, 0, {'Approx.';'d_n_o_r_m'}, 'HorizontalAlignment', 'center');
            
            line([min(arrayXRangeInterp) max(arrayXRangeInterp)], [max(arrayYRangeInterp) max(arrayYRangeInterp)], [0 0], 'LineWidth', numOverlayBorderWidth, 'Color', [1 0.5 0], 'LineStyle', '-');
            text(min(arrayXRangeInterp),min(arrayYRangeInterp),max(imageSubAreaInterp(:))*1.5,'C','FontSize',numSubFigLabelFontSize,'Color','k', 'FontWeight', 'bold');
            
        elseif structOtherImageOne(iOutputProtein).dimensionality == 3,
            %volume plot
            hold on;
            
            if isdeployed,
                %don't use alpha data or lighting, because this can crash 
                % the virtual box (due to OpenGL rendering issues)
                plotPatch1 = patch(structOtherImageOne(iOutputProtein).fv1, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none' );

                plotPatch2 = patch(structOtherImageOne(iOutputProtein).fv2, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none' );

                plotPatch3 = patch(structOtherImageOne(iOutputProtein).fv3, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none' );
            else
                %plot using alpha transparency
                plotPatch1 = patch(structOtherImageOne(iOutputProtein).fv1, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none', ...
                    'FaceAlpha', 0.05, ...
                    'FaceLighting','phong', ...
                    'AmbientStrength',0.1);

                plotPatch2 = patch(structOtherImageOne(iOutputProtein).fv2, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none', ...
                    'FaceAlpha', 0.5, ...
                    'FaceLighting','phong', ...
                    'AmbientStrength',0.5);

                plotPatch3 = patch(structOtherImageOne(iOutputProtein).fv3, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none', ...
                    'FaceAlpha', 0.7, ...
                    'FaceLighting','phong', ...
                    'AmbientStrength', 0.7);
                %and use lighting
                camlight;
            end
            
            line( [0 arrayOtherImageCoOrds{iOutputProtein,1}(2)-arrayOtherImageCoOrds{iOutputProtein,1}(1)], ...
                  [arrayOtherImageCoOrds{iOutputProtein,1}(4)-arrayOtherImageCoOrds{iOutputProtein,1}(3) arrayOtherImageCoOrds{iOutputProtein,1}(4)-arrayOtherImageCoOrds{iOutputProtein,1}(3)], ...
                  [0 0], 'LineWidth', numOverlayBorderWidth, 'Color', [1 0.5 0], 'LineStyle', '-');
            
            box on;
            axis tight;
            axis ij;
            view(arrayOtherImageView{iOutputProtein,1});
            axis equal;
            set(gca, 'Projection', 'Perspective');
            set(gca, 'YTick', [], 'YTickLabel', [], 'ZTick', [], 'ZTickLabel', []);
            set(gca, 'XTick', [0 arrayOtherImageCoOrds{iOutputProtein,1}(2)-arrayOtherImageCoOrds{iOutputProtein,1}(1)], 'XTickLabel', [{['~' num2str(arrayOutSampleAnalysis(1).NormDist, '%3.2f')]}; {['~' num2str(arrayOutSampleAnalysis(2).NormDist, '%3.2f')]}]);
            xlabel('Approx. d_n_o_r_m');
            set(gca, 'DataAspectRatio', [1 1 0.6]);
            hold off;
            text(0,0,(arrayOtherImageCoOrds{iOutputProtein,1}(6)-arrayOtherImageCoOrds{iOutputProtein,1}(5))*1.5,'C','FontSize',numSubFigLabelFontSize,'Color','k', 'FontWeight', 'bold');
            
        end
        
        %plot the second 'other image'
        handleOtherImage = subplot('Position', arrayTwoLocOtherImagePositionDoubleTwo);
        if structOtherImageTwo(iOutputProtein).dimensionality == 2,
            %surface plot
            arrayXRange = structOtherImageOne(iOutputProtein).arrayXCoOrds;
            arrayYRange = structOtherImageOne(iOutputProtein).arrayYCoOrds;
            arrayXRangeInterp = arrayXRange(1):(1/double(arrayNumSurfaceInterp(iOutputProtein,1))):arrayXRange(end);
            arrayYRangeInterp = arrayYRange(1):(1/double(arrayNumSurfaceInterp(iOutputProtein,1))):arrayYRange(end);
            imageSubArea = structOtherImageOne(iOutputProtein).arrayPixelInt;
            imageSubAreaInterp = interp2(arrayXRange, arrayYRange, ...
                                         double(imageSubArea(:,:)), ...
                                         arrayXRangeInterp, arrayYRangeInterp', ...
                                         'cubic');
            arraySubAreaInterpBelowZeroObsIndices = find(imageSubAreaInterp < 0.01);
            imageSubAreaInterp(arraySubAreaInterpBelowZeroObsIndices) = 0.01;
            if isdeployed,
                %don't use alpha data or lighting, because this can crash 
                % the virtual box (due to OpenGL rendering issues)
                surf(arrayXRangeInterp, arrayYRangeInterp, imageSubAreaInterp, ...
                    'FaceColor','interp',  'EdgeColor','none' );
            else
                %plot using alpha transparency
                arrayAlphaData = log(imageSubAreaInterp);

                surf(arrayXRangeInterp, arrayYRangeInterp, imageSubAreaInterp, ...
                    'FaceColor','interp',  'EdgeColor','none',  'FaceAlpha','interp',...
                    'AlphaDataMapping','scaled',  'AlphaData',arrayAlphaData, ...
                    'FaceLighting','phong', 'AmbientStrength',0.8);
                %and use lighting
                camlight;
            end
            axis tight;
            axis ij;
            axis equal;
            colormap(jet);
            set(gca, 'DataAspectRatio',[1 1 4/double(arrayNumSurfaceInterp(iOutputProtein,1))]);
            view(arrayOtherImageView{iOutputProtein,1}(1), arrayOtherImageView{iOutputProtein,1}(2));
            set(gca, 'YTick', [], 'YTickLabel', [], 'ZTick', [], 'ZTickLabel', []);
            set(gca, 'XTick', [arrayXRangeInterp(1) arrayXRangeInterp(end)], 'XTickLabel', [{['~' num2str(arrayOutSampleAnalysis(1).NormDist, '%3.2f')]}; {['~' num2str(arrayOutSampleAnalysis(2).NormDist, '%3.2f')]}]);
            zlabel('Pix. Int.');
            text((arrayOtherImageCoOrds{iOutputProtein,1}(2)-arrayOtherImageCoOrds{iOutputProtein,1}(1))/2, (arrayOtherImageCoOrds{iOutputProtein,1}(4)-arrayOtherImageCoOrds{iOutputProtein,1}(3))*1.35, 0, {'Approx.';'d_n_o_r_m'}, 'HorizontalAlignment', 'center');
            line([min(arrayXRangeInterp) max(arrayXRangeInterp)], [max(arrayYRangeInterp) max(arrayYRangeInterp)], [0 0], 'LineWidth', numOverlayBorderWidth, 'Color', [1 0.5 0], 'LineStyle', '-');
            text(min(arrayXRangeInterp),min(arrayYRangeInterp),max(imageSubAreaInterp(:))*1.5,'D','FontSize',numSubFigLabelFontSize,'Color','k', 'FontWeight', 'bold');
            
        elseif structOtherImageTwo(iOutputProtein).dimensionality == 3,
            %volume plot
            hold on;
            if isdeployed,
                %don't use alpha data or lighting, because this can crash 
                % the virtual box (due to OpenGL rendering issues)
                plotPatch1 = patch(structOtherImageTwo(iOutputProtein).fv1, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none' );

                plotPatch2 = patch(structOtherImageTwo(iOutputProtein).fv2, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none' );

                plotPatch3 = patch(structOtherImageTwo(iOutputProtein).fv3, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none' );
            else
                
                plotPatch1 = patch(structOtherImageTwo(iOutputProtein).fv1, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none', ...
                    'FaceAlpha', 0.05, ...
                    'FaceLighting','phong', ...
                    'AmbientStrength',0.1);

                plotPatch2 = patch(structOtherImageTwo(iOutputProtein).fv2, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none', ...
                    'FaceAlpha', 0.5, ...
                    'FaceLighting','phong', ...
                    'AmbientStrength',0.5);

                plotPatch3 = patch(structOtherImageTwo(iOutputProtein).fv3, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none', ...
                    'FaceAlpha', 0.7, ...
                    'FaceLighting','phong', ...
                    'AmbientStrength', 0.7);

                camlight;
            end
            
            line( [0 arrayOtherImageCoOrds{iOutputProtein,2}(2)-arrayOtherImageCoOrds{iOutputProtein,2}(1)], ...
                  [arrayOtherImageCoOrds{iOutputProtein,2}(4)-arrayOtherImageCoOrds{iOutputProtein,2}(3) arrayOtherImageCoOrds{iOutputProtein,2}(4)-arrayOtherImageCoOrds{iOutputProtein,2}(3)], ...
                  [0 0], 'LineWidth', numOverlayBorderWidth, 'Color', [1 0.5 0], 'LineStyle', '-');
            
            box on;
            axis tight;
            axis ij;
            view(arrayOtherImageView{iOutputProtein,2});
            axis equal;
            set(gca, 'Projection', 'Perspective');
            set(gca, 'YTick', [], 'YTickLabel', [], 'ZTick', [], 'ZTickLabel', []);
            set(gca, 'XTick', [0 arrayOtherImageCoOrds{iOutputProtein,2}(2)-arrayOtherImageCoOrds{iOutputProtein,2}(1)], 'XTickLabel', [{['~' num2str(arrayOutSampleAnalysis(1).NormDist, '%3.2f')]}; {['~' num2str(arrayOutSampleAnalysis(2).NormDist, '%3.2f')]}]);
            text((arrayOtherImageCoOrds{iOutputProtein,2}(2)-arrayOtherImageCoOrds{iOutputProtein,2}(1))/2, (arrayOtherImageCoOrds{iOutputProtein,2}(4)-arrayOtherImageCoOrds{iOutputProtein,2}(3))*1.35, 0, {'Approx.';'d_n_o_r_m'}, 'HorizontalAlignment', 'center');
            set(gca, 'DataAspectRatio', [1 1 0.6]);
            hold off;
            text(0,0,(arrayOtherImageCoOrds{iOutputProtein,2}(6)-arrayOtherImageCoOrds{iOutputProtein,2}(5))*1.8,'D','FontSize',numSubFigLabelFontSize,'Color','k', 'FontWeight', 'bold');
            
        end
        
    elseif numLocalisationsForProtein == 3,
        
        %plot the first 'other image'
        handleOtherImage = subplot('Position', arrayThreeLocOtherImagePositionOne);
        if structOtherImageOne(iOutputProtein).dimensionality == 2,
            %surface plot
            arrayXRange = structOtherImageOne(iOutputProtein).arrayXCoOrds;
            arrayYRange = structOtherImageOne(iOutputProtein).arrayYCoOrds;
            arrayXRangeInterp = arrayXRange(1):(1/double(arrayNumSurfaceInterp(iOutputProtein,1))):arrayXRange(end);
            arrayYRangeInterp = arrayYRange(1):(1/double(arrayNumSurfaceInterp(iOutputProtein,1))):arrayYRange(end);
            imageSubArea = structOtherImageOne(iOutputProtein).arrayPixelInt;
            imageSubAreaInterp = interp2(arrayXRange, arrayYRange, ...
                                         double(imageSubArea(:,:)), ...
                                         arrayXRangeInterp, arrayYRangeInterp', ...
                                         'cubic');
            arraySubAreaInterpBelowZeroObsIndices = find(imageSubAreaInterp < 0.01);
            imageSubAreaInterp(arraySubAreaInterpBelowZeroObsIndices) = 0.01;
            if isdeployed,
                %don't use alpha data or lighting, because this can crash 
                % the virtual box (due to OpenGL rendering issues)
                surf(arrayXRangeInterp, arrayYRangeInterp, imageSubAreaInterp, ...
                    'FaceColor','interp',  'EdgeColor','none');
            else

                arrayAlphaData = log(imageSubAreaInterp);

                surf(arrayXRangeInterp, arrayYRangeInterp, imageSubAreaInterp, ...
                    'FaceColor','interp',  'EdgeColor','none',  'FaceAlpha','interp',...
                    'AlphaDataMapping','scaled',  'AlphaData',arrayAlphaData, ...
                    'FaceLighting','phong', 'AmbientStrength',0.8);
                
                camlight;
            end
            axis tight;
            axis ij;
            axis equal;
            colormap(jet);
            set(gca, 'DataAspectRatio',[1 1 4/double(arrayNumSurfaceInterp(iOutputProtein,1))]);
            view(arrayOtherImageView{iOutputProtein,1}(1), arrayOtherImageView{iOutputProtein,1}(2));
            set(gca, 'YTick', [], 'YTickLabel', [], 'ZTick', [], 'ZTickLabel', []);
            set(gca, 'XTick', [arrayXRangeInterp(1) arrayXRangeInterp(end)], 'XTickLabel', [{['~' num2str(arrayOutSampleAnalysis(1).NormDist, '%3.2f')]}; {['~' num2str(arrayOutSampleAnalysis(2).NormDist, '%3.2f')]}]);
            text((arrayXRange(end) + arrayXRange(1))/2, arrayYRange(end)*1.35, 0, {'Approx.';'d_n_o_r_m'}, 'HorizontalAlignment', 'center');
            zlabel('Pix. Int.');
            line([min(arrayXRangeInterp) max(arrayXRangeInterp)], [max(arrayYRangeInterp) max(arrayYRangeInterp)], [0 0], 'LineWidth', numOverlayBorderWidth, 'Color', [1 0.5 0], 'LineStyle', '-');
            text(min(arrayXRangeInterp),min(arrayYRangeInterp),max(imageSubAreaInterp(:))*1.5,'C','FontSize',numSubFigLabelFontSize,'Color','k', 'FontWeight', 'bold');
            
        elseif structOtherImageOne(iOutputProtein).dimensionality == 3,
            %volume plot
            hold on;
            if isdeployed,
                %don't use alpha data or lighting, because this can crash 
                % the virtual box (due to OpenGL rendering issues)
                plotPatch1 = patch(structOtherImageOne(iOutputProtein).fv1, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none' );

                plotPatch2 = patch(structOtherImageOne(iOutputProtein).fv2, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none' );

                plotPatch3 = patch(structOtherImageOne(iOutputProtein).fv3, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none' );
            else
                plotPatch1 = patch(structOtherImageOne(iOutputProtein).fv1, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none', ...
                    'FaceAlpha', 0.05, ...
                    'FaceLighting','phong', ...
                    'AmbientStrength',0.1);

                plotPatch2 = patch(structOtherImageOne(iOutputProtein).fv2, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none', ...
                    'FaceAlpha', 0.5, ...
                    'FaceLighting','phong', ...
                    'AmbientStrength',0.5);

                plotPatch3 = patch(structOtherImageOne(iOutputProtein).fv3, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none', ...
                    'FaceAlpha', 0.7, ...
                    'FaceLighting','phong', ...
                    'AmbientStrength', 0.7);
                
                camlight;
            end
            
            line( [0 arrayOtherImageCoOrds{iOutputProtein,1}(2)-arrayOtherImageCoOrds{iOutputProtein,1}(1)], ...
                  [arrayOtherImageCoOrds{iOutputProtein,1}(4)-arrayOtherImageCoOrds{iOutputProtein,1}(3) arrayOtherImageCoOrds{iOutputProtein,1}(4)-arrayOtherImageCoOrds{iOutputProtein,1}(3)], ...
                  [0 0], 'LineWidth', numOverlayBorderWidth, 'Color', [1 0.5 0], 'LineStyle', '-');
            
            box on;
            axis tight;
            axis ij;
            view(arrayOtherImageView{iOutputProtein,1});
            axis equal;
            set(gca, 'Projection', 'Perspective');
            set(gca, 'YTick', [], 'YTickLabel', [], 'ZTick', [], 'ZTickLabel', []);
            set(gca, 'XTick', [0 arrayOtherImageCoOrds{iOutputProtein,1}(2)-arrayOtherImageCoOrds{iOutputProtein,1}(1)], 'XTickLabel', [{['~' num2str(arrayOutSampleAnalysis(1).NormDist, '%3.2f')]}; {['~' num2str(arrayOutSampleAnalysis(2).NormDist, '%3.2f')]}]);
            xlabel('Approx. d_n_o_r_m');
            set(gca, 'DataAspectRatio', [1 1 0.6]);
            hold off;
            text(0,0,(arrayOtherImageCoOrds{iOutputProtein,1}(6)-arrayOtherImageCoOrds{iOutputProtein,1}(5))*1.5,'C','FontSize',numSubFigLabelFontSize,'Color','k', 'FontWeight', 'bold');
            
        end
        
        
        %plot the first 'other image'
        handleOtherImage = subplot('Position', arrayThreeLocOtherImagePositionTwo);
        if structOtherImageTwo(iOutputProtein).dimensionality == 2,
            %surface plot
            arrayXRange = structOtherImageTwo(iOutputProtein).arrayXCoOrds;
            arrayYRange = structOtherImageTwo(iOutputProtein).arrayYCoOrds;
            arrayXRangeInterp = arrayXRange(1):(1/double(arrayNumSurfaceInterp(iOutputProtein,1))):arrayXRange(end);
            arrayYRangeInterp = arrayYRange(1):(1/double(arrayNumSurfaceInterp(iOutputProtein,1))):arrayYRange(end);
            imageSubArea = structOtherImageTwo(iOutputProtein).arrayPixelInt;
            imageSubAreaInterp = interp2(arrayXRange, arrayYRange, ...
                                         double(imageSubArea(:,:)), ...
                                         arrayXRangeInterp, arrayYRangeInterp', ...
                                         'cubic');
            arraySubAreaInterpBelowZeroObsIndices = find(imageSubAreaInterp < 0.01);
            imageSubAreaInterp(arraySubAreaInterpBelowZeroObsIndices) = 0.01;

            if isdeployed,
                %don't use alpha data or lighting, because this can crash 
                % the virtual box (due to OpenGL rendering issues)
                surf(arrayXRangeInterp, arrayYRangeInterp, imageSubAreaInterp, ...
                    'FaceColor','interp',  'EdgeColor','none' );
            else
                arrayAlphaData = log(imageSubAreaInterp);
                surf(arrayXRangeInterp, arrayYRangeInterp, imageSubAreaInterp, ...
                    'FaceColor','interp',  'EdgeColor','none',  'FaceAlpha','interp',...
                    'AlphaDataMapping','scaled',  'AlphaData',arrayAlphaData, ...
                    'FaceLighting','phong', 'AmbientStrength',0.8);
                camlight;
            end
            axis tight;
            axis ij;
            axis equal;
            
            colormap(jet);
            set(gca, 'DataAspectRatio',[1 1 4/double(arrayNumSurfaceInterp(iOutputProtein,1))]);
            view(arrayOtherImageView{iOutputProtein,1}(1), arrayOtherImageView{iOutputProtein,1}(2));
            set(gca, 'YTick', [], 'YTickLabel', [], 'ZTick', [], 'ZTickLabel', []);
            set(gca, 'XTick', [arrayXRangeInterp(1) arrayXRangeInterp(end)], 'XTickLabel', [{['~' num2str(arrayOutSampleAnalysis(1).NormDist, '%3.2f')]}; {['~' num2str(arrayOutSampleAnalysis(2).NormDist, '%3.2f')]}]);
            xlabel('Approx. d_n_o_r_m');
            zlabel('Pix. Int.');
            line([min(arrayXRangeInterp) max(arrayXRangeInterp)], [max(arrayYRangeInterp) max(arrayYRangeInterp)], [0 0], 'LineWidth', numOverlayBorderWidth, 'Color', [1 0.5 0], 'LineStyle', '-');
            text(min(arrayXRangeInterp),min(arrayYRangeInterp),max(imageSubAreaInterp(:))*1.5,'D','FontSize',numSubFigLabelFontSize,'Color','k', 'FontWeight', 'bold');
            
        elseif structOtherImageTwo(iOutputProtein).dimensionality == 3,
            %volume plot
            hold on;
            if isdeployed,
                %don't use alpha data or lighting, because this can crash 
                % the virtual box (due to OpenGL rendering issues)
                plotPatch1 = patch(structOtherImageTwo(iOutputProtein).fv1, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none' );

                plotPatch2 = patch(structOtherImageTwo(iOutputProtein).fv2, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none' );

                plotPatch3 = patch(structOtherImageTwo(iOutputProtein).fv3, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none' );
            else
                
                plotPatch1 = patch(structOtherImageTwo(iOutputProtein).fv1, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none', ...
                    'FaceAlpha', 0.05, ...
                    'FaceLighting','phong', ...
                    'AmbientStrength',0.1);

                plotPatch2 = patch(structOtherImageTwo(iOutputProtein).fv2, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none', ...
                    'FaceAlpha', 0.5, ...
                    'FaceLighting','phong', ...
                    'AmbientStrength',0.5);

                plotPatch3 = patch(structOtherImageTwo(iOutputProtein).fv3, ... 
                    'FaceColor','g', ... 
                    'edgecolor', 'none', ...
                    'FaceAlpha', 0.7, ...
                    'FaceLighting','phong', ...
                    'AmbientStrength', 0.7);
                
                camlight;
            end
            line( [0 arrayOtherImageCoOrds{iOutputProtein,2}(2)-arrayOtherImageCoOrds{iOutputProtein,2}(1)], ...
                  [arrayOtherImageCoOrds{iOutputProtein,2}(4)-arrayOtherImageCoOrds{iOutputProtein,2}(3) arrayOtherImageCoOrds{iOutputProtein,2}(4)-arrayOtherImageCoOrds{iOutputProtein,2}(3)], ...
                  [0 0], 'LineWidth', numOverlayBorderWidth, 'Color', [1 0.5 0], 'LineStyle', '-');
            
            box on;
            axis tight;
            axis ij;
            view(arrayOtherImageView{iOutputProtein,2});
            axis equal;
            set(gca, 'Projection', 'Perspective');
            set(gca, 'YTick', [], 'YTickLabel', [], 'ZTick', [], 'ZTickLabel', []);
            set(gca, 'XTick', [0 arrayOtherImageCoOrds{iOutputProtein,1}(2)-arrayOtherImageCoOrds{iOutputProtein,1}(1)], 'XTickLabel', [{['~' num2str(arrayOutSampleAnalysis(1).NormDist, '%3.2f')]}; {['~' num2str(arrayOutSampleAnalysis(2).NormDist, '%3.2f')]}]);
            xlabel('Approx. d_n_o_r_m');
            set(gca, 'DataAspectRatio', [1 1 0.6]);
            hold off;
            text(0,0,(arrayOtherImageCoOrds{iOutputProtein,1}(6)-arrayOtherImageCoOrds{iOutputProtein,1}(5))*1.5,'D','FontSize',numSubFigLabelFontSize,'Color','k', 'FontWeight', 'bold');
            
        end
        
    end

    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% Plot the lowess-smoothed curves
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%     for iLocalisation = 1:numLocalisationsForProtein,
                
    %plot the membrane localisation first (if it exists)
    if ~(arrayGroupedNodes(iOutputProtein,3) == 0),

        numNode = arrayGroupedNodes(iOutputProtein,3);

        %load the data cloud for the representative image
        stringRepImageMembDataFolder = [ stringRepImageDataFolder 'M' strFolderSep ];

        structRepImageSampleData = loadSampAnalysis(stringRepImageMembDataFolder, '.mat');

        numObjects = length(structRepImageSampleData);
        numSamplesPerObject = size(structRepImageSampleData(1).NormDist,1);
        numPixelsPerSample = size(structRepImageSampleData(1).SigInt,2);
        arraySampleValues = zeros(numObjects*numSamplesPerObject, numPixelsPerSample, 'uint8');
        arraySamplePositions = zeros(numObjects*numSamplesPerObject,1, 'double');
        for iObject = 1:numObjects,
            for iSample = 1:numSamplesPerObject,
                arraySamplePositions(  (iObject-1)*numSamplesPerObject + iSample  ) = structRepImageSampleData(iObject).NormDist(iSample);
                arraySampleValues(  (iObject-1)*numSamplesPerObject + iSample, :  ) = structRepImageSampleData(iObject).SigInt(iSample,:);
            end 
        end
        arraySamplePositions = rescaleNormDist(arraySamplePositions, numSpatialBins, arrayNormDistRatio);
        numRepImageSmpMean = mean(double(arraySampleValues(:)));
        numRepImageStDev = std(double(arraySampleValues(:)));


        handlePlotOne = subplot('Position', arrayPlot1Position);


        %overlay the lowess smoothed curves and the data cloud for the representative image
        hold on;
        plot(arraySamplePositions, (double(arraySampleValues) - numRepImageSmpMean)/numRepImageStDev, '.', 'Color', arrayDataCloudColor(arrayPatientToDisplay(iOutputProtein),:), 'MarkerSize', numPlotDataCloudMarkerSize);


        for iPatient = 1:numPatients,

            stringMembDataFolder = [ strProcDataFolder arrayOutputProtFolders{iOutputProtein} strFolderSep 'Pat_' num2str(iPatient) strFolderSep 'M' strFolderSep ];
            structSampleData = loadSampAnalysis(stringMembDataFolder, '.mat');
            arrayLoessData = loadLoessCurve(stringMembDataFolder, '.mat');

            numDataObjects = length(structSampleData);
            numSamplesPerObject = length(structSampleData(1).SigInt);
            arraySampleData  = zeros(numDataObjects*numSamplesPerObject,1,'uint8');
            for iObj = 1:numDataObjects,
                numBaseIndex = (iObj-1)*numSamplesPerObject;

                arraySampleData((numBaseIndex+1):(numBaseIndex+iSample)) = structSampleData(iObj).SigInt(:);

            end

            numDataMean = mean(double(arraySampleData));
            numDataStDev = std(double(arraySampleData));

            plot(arrayLoessData(1,:), (arrayLoessData(2,:) - numDataMean)/numDataStDev, '.', 'Color', arrayPatientColors{iPatient}, 'MarkerSize', numPlotLowessMarkerSize);
            arrayUpperLoessCI = [arrayLoessCIBounds{1, numNode, 1}(:,1); arrayLoessCIBounds{1, numNode, 2}(:,1); arrayLoessCIBounds{1, numNode, 3}(:,1)];
            arrayLowerLoessCI = [arrayLoessCIBounds{1, numNode, 1}(:,2); arrayLoessCIBounds{1, numNode, 2}(:,2); arrayLoessCIBounds{1, numNode, 3}(:,2)];
            arrayLoessCIXPosCombined = [arrayLoessCIXPos{1, numNode, 1}(:); arrayLoessCIXPos{1, numNode, 2}(:); arrayLoessCIXPos{1, numNode, 3}(:)];

            arrayLoessX = arrayLoessData(1,:);
            arrayLoessVals = arrayLoessData(2,:);


            for iTissueLayer = 1:numTissueLayers,

                numTissueLayerMinX = double(arrayDivisionIndices(iTissueLayer))-1;
                numTissueLayerMaxX = double(arrayDivisionIndices(iTissueLayer+1))-1;

                if iTissueLayer < numTissueLayers,
                    arrayIndexInTissueLayer = find( (arrayLoessX >= numTissueLayerMinX) & (arrayLoessX < numTissueLayerMaxX) );
                    arrayCIIndexInTissueLayer = find( (arrayLoessCIXPosCombined >= numTissueLayerMinX) & (arrayLoessCIXPosCombined < numTissueLayerMaxX) );
                else
                    arrayIndexInTissueLayer = find( (arrayLoessX >= numTissueLayerMinX) & (arrayLoessX <= numTissueLayerMaxX) );
                    arrayCIIndexInTissueLayer = find( (arrayLoessCIXPosCombined >= numTissueLayerMinX) & (arrayLoessCIXPosCombined <= numTissueLayerMaxX) );
                end
                arrayXinTissueLayer = arrayLoessX(arrayIndexInTissueLayer);
                arrayCIXinTissueLayer = arrayLoessCIXPosCombined(arrayCIIndexInTissueLayer);

                numPlotMinX = max([arrayXinTissueLayer(1) arrayCIXinTissueLayer(1)]);
                numPlotMaxX = min([arrayXinTissueLayer(end) arrayCIXinTissueLayer(end)]);


                arrayXToPlot = linspace(numPlotMinX, numPlotMaxX,40);

                arrayLoessToPlot = interp1(arrayXinTissueLayer, arrayLoessVals(arrayIndexInTissueLayer), arrayXToPlot);
                arrayUpperLoessToPlot = interp1(arrayCIXinTissueLayer, arrayUpperLoessCI(arrayCIIndexInTissueLayer), arrayXToPlot);
                arrayLowerLoessToPlot = interp1(arrayCIXinTissueLayer, arrayLowerLoessCI(arrayCIIndexInTissueLayer), arrayXToPlot);


                plot(arrayXToPlot, (arrayLoessToPlot - numDataMean + arrayUpperLoessToPlot)/numDataStDev, '--', 'Color', arrayPatientColors{iPatient},'LineWidth',(numPlotLowessMarkerSize/8));
                plot(arrayXToPlot, (arrayLoessToPlot - numDataMean + arrayLowerLoessToPlot)/numDataStDev, '--', 'Color', arrayPatientColors{iPatient},'LineWidth',(numPlotLowessMarkerSize/8));

            end
        end

        hold off;


        %extract the max and minimum values for adjusting axes
        arrayTempMax = (double(max(arraySampleValues(:))) - numRepImageSmpMean)/numRepImageStDev;
        arrayTempMin = (double(min(arraySampleValues(:))) - numRepImageSmpMean)/numRepImageStDev;     
        if arrayTempMax > 3 || isnan(arrayTempMax),
            arrayTempMax = 3;
        end
        if arrayTempMin < -3 || isnan(arrayTempMin),
            arrayTempMin = -3;
        end

        %set the axes
        axis([0 numSpatialBins, arrayTempMin*1.2 arrayTempMax*1.2]);
        text(  numSpatialBins*0.75, arrayTempMax*0.85, {'\bfPlasma';'\bfMembrane'},  'HorizontalAlignment', 'center',  'FontSize', numPlotFontSizeTitle  );
        set(gca, 'XTick', arrayLayerBoundaries, 'XTickLabel', [0 1 2 3]);
        set(gca, 'YTick', [ -4, -3, -2, -1, 0, 1, 2, 3, 4 ], 'YTickLabel', []);
        for iTick = ceil(arrayTempMin):floor(arrayTempMax),
            if (iTick >= 2),
                text(-1, double(iTick), ['$$\bar{\mu}$$+' num2str(iTick) '$$\bar{\sigma}$$'], 'Interpreter', 'latex', 'HorizontalAlignment', 'right');
            elseif (iTick <= -2),
                text(-1, double(iTick), ['$$\bar{\mu}$$-' num2str(iTick) '$$\bar{\sigma}$$'], 'Interpreter', 'latex', 'HorizontalAlignment', 'right');
            elseif iTick == 0,
                text(-1, double(iTick), '$$\bar{\mu}$$', 'Interpreter', 'latex', 'HorizontalAlignment', 'right');
            elseif iTick == 1
                text(-1, double(iTick), ['$$\bar{\mu}$$+$$\bar{\sigma}$$'], 'Interpreter', 'latex', 'HorizontalAlignment', 'right');
            elseif iTick == -1
                text(-1, double(iTick), '$$\bar{\mu}$$-$$-\bar{\sigma}$$', 'Interpreter', 'latex', 'HorizontalAlignment', 'right');
            end
        end
        set(gca, 'FontSize', numPlotFontSize);
        xlabel('d_n_o_r_m','FontSize',numPlotFontSizeTitle, 'Position', [14 arrayTempMin*1.35 1]);
        ylabel('I_s_,_n_o_r_m','FontSize',numPlotFontSizeTitle, 'Position', [-4.5 (arrayTempMax+arrayTempMin)/2 1]);
        if (  ( numLocalisationsForProtein == 1 ) ),
            text(-7.5,arrayTempMax*0.95,'D','FontSize',numSubFigLabelFontSize,'Color','k', 'FontWeight', 'bold');
        elseif (numLocalisationsForProtein == 2) && (length(arrayOtherImageTypes{iOutputProtein}) == 1),
            text(-7.5,arrayTempMax*0.95,'D','FontSize',numSubFigLabelFontSize,'Color','k', 'FontWeight', 'bold');
        elseif (numLocalisationsForProtein == 2) && (length(arrayOtherImageTypes{iOutputProtein}) == 2),
            text(-7.5,arrayTempMax*0.95,'E','FontSize',numSubFigLabelFontSize,'Color','k', 'FontWeight', 'bold');
        elseif numLocalisationsForProtein == 3,
            text(-7.5,arrayTempMax*0.95,'E','FontSize',numSubFigLabelFontSize,'Color','k', 'FontWeight', 'bold');
        end

        %the cytoplasmic localisation is always sampled, plot it second
        if ~(arrayGroupedNodes(iOutputProtein,1) == 0),
            numNode = arrayGroupedNodes(iOutputProtein,1);

            %load the data cloud for the representative image
            stringRepImageCytoDataFolder = [ stringRepImageDataFolder 'C' strFolderSep ];

            structRepImageSampleData = loadSampAnalysis(stringRepImageCytoDataFolder, '.mat');

            numObjects = length(structRepImageSampleData);
            numSamplesPerObject = size(structRepImageSampleData(1).NormDist,1);
            numPixelsPerSample = size(structRepImageSampleData(1).SigInt,2);
            arraySampleValues = zeros(numObjects*numSamplesPerObject, numPixelsPerSample, 'uint8');
            arraySamplePositions = zeros(numObjects*numSamplesPerObject,1, 'double');
            for iObject = 1:numObjects,
                for iSample = 1:numSamplesPerObject,
                    arraySamplePositions(  (iObject-1)*numSamplesPerObject + iSample  ) = structRepImageSampleData(iObject).NormDist(iSample);
                    arraySampleValues(  (iObject-1)*numSamplesPerObject + iSample, :  ) = structRepImageSampleData(iObject).SigInt(iSample,:);
                end 
            end
            arraySamplePositions = rescaleNormDist(arraySamplePositions, numSpatialBins, arrayNormDistRatio);
            numRepImageSmpMean = mean(double(arraySampleValues(:)));
            numRepImageStDev = std(double(arraySampleValues(:)));


            handlePlotTwo = subplot('Position', arrayPlot2Position);

            %overlay the lowess smoothed curves and the data cloud for
            %the representative image
            hold on;
            plot(arraySamplePositions, (double(arraySampleValues) - numRepImageSmpMean)/numRepImageStDev, '.', 'Color', arrayDataCloudColor(arrayPatientToDisplay(iOutputProtein),:), 'MarkerSize', numPlotDataCloudMarkerSize);


            for iPatient = 1:numPatients,

                stringCytoDataFolder = [ strProcDataFolder arrayOutputProtFolders{iOutputProtein} strFolderSep 'Pat_' num2str(iPatient) strFolderSep 'C' strFolderSep ];
                structSampleData = loadSampAnalysis(stringCytoDataFolder, '.mat');
                arrayLoessData = loadLoessCurve(stringCytoDataFolder, '.mat');

                numDataObjects = length(structSampleData);
                numSamplesPerObject = length(structSampleData(1).SigInt);
                arraySampleData  = zeros(numDataObjects*numSamplesPerObject,1,'uint8');
                for iObj = 1:numDataObjects,
                    numBaseIndex = (iObj-1)*numSamplesPerObject;

                    arraySampleData((numBaseIndex+1):(numBaseIndex+iSample)) = structSampleData(iObj).SigInt(:);

                end

                numDataMean = mean(double(arraySampleData));
                numDataStDev = std(double(arraySampleData));



                plot(arrayLoessData(1,:), (arrayLoessData(2,:) - numDataMean)/numDataStDev, '.', 'Color', arrayPatientColors{iPatient}, 'MarkerSize', numPlotLowessMarkerSize);
                arrayUpperLoessCI = [arrayLoessCIBounds{1, numNode, 1}(:,1); arrayLoessCIBounds{1, numNode, 2}(:,1); arrayLoessCIBounds{1, numNode, 3}(:,1)];
                arrayLowerLoessCI = [arrayLoessCIBounds{1, numNode, 1}(:,2); arrayLoessCIBounds{1, numNode, 2}(:,2); arrayLoessCIBounds{1, numNode, 3}(:,2)];

                arrayLoessCIXPosCombined = [arrayLoessCIXPos{1, numNode, 1}(:); arrayLoessCIXPos{1, numNode, 2}(:); arrayLoessCIXPos{1, numNode, 3}(:)];

                arrayLoessX = arrayLoessData(1,:);
                arrayLoessVals = arrayLoessData(2,:);


                for iTissueLayer = 1:numTissueLayers,

                    numTissueLayerMinX = double(arrayDivisionIndices(iTissueLayer))-1;
                    numTissueLayerMaxX = double(arrayDivisionIndices(iTissueLayer+1))-1;

                    if iTissueLayer < numTissueLayers,
                        arrayIndexInTissueLayer = find( (arrayLoessX >= numTissueLayerMinX) & (arrayLoessX < numTissueLayerMaxX) );
                        arrayCIIndexInTissueLayer = find( (arrayLoessCIXPosCombined >= numTissueLayerMinX) & (arrayLoessCIXPosCombined < numTissueLayerMaxX) );
                    else
                        arrayIndexInTissueLayer = find( (arrayLoessX >= numTissueLayerMinX) & (arrayLoessX <= numTissueLayerMaxX) );
                        arrayCIIndexInTissueLayer = find( (arrayLoessCIXPosCombined >= numTissueLayerMinX) & (arrayLoessCIXPosCombined <= numTissueLayerMaxX) );
                    end
                    arrayXinTissueLayer = arrayLoessX(arrayIndexInTissueLayer);
                    arrayCIXinTissueLayer = arrayLoessCIXPosCombined(arrayCIIndexInTissueLayer);

                    numPlotMinX = max([arrayXinTissueLayer(1) arrayCIXinTissueLayer(1)]);
                    numPlotMaxX = min([arrayXinTissueLayer(end) arrayCIXinTissueLayer(end)]);


                    arrayXToPlot = linspace(numPlotMinX, numPlotMaxX,40);

                    arrayLoessToPlot = interp1(arrayXinTissueLayer, arrayLoessVals(arrayIndexInTissueLayer), arrayXToPlot);
                    arrayUpperLoessToPlot = interp1(arrayCIXinTissueLayer, arrayUpperLoessCI(arrayCIIndexInTissueLayer), arrayXToPlot);
                    arrayLowerLoessToPlot = interp1(arrayCIXinTissueLayer, arrayLowerLoessCI(arrayCIIndexInTissueLayer), arrayXToPlot);


                    plot(arrayXToPlot, (arrayLoessToPlot - numDataMean + arrayUpperLoessToPlot)/numDataStDev, '--', 'Color', arrayPatientColors{iPatient},'LineWidth',(numPlotLowessMarkerSize/8));
                    plot(arrayXToPlot, (arrayLoessToPlot - numDataMean + arrayLowerLoessToPlot)/numDataStDev, '--', 'Color', arrayPatientColors{iPatient},'LineWidth',(numPlotLowessMarkerSize/8));

                end
            end

            hold off;



            %extract the max and minimum values for adjusting axes
            arrayTempMax = (double(max(arraySampleValues(:))) - numRepImageSmpMean)/numRepImageStDev;
            arrayTempMin = (double(min(arraySampleValues(:))) - numRepImageSmpMean)/numRepImageStDev;     
            if arrayTempMax > 3 || isnan(arrayTempMax),
                arrayTempMax = 3;
            end
            if arrayTempMin < -3 || isnan(arrayTempMin),
                arrayTempMin = -3;
            end

            %set the axes
            axis([0 numSpatialBins, arrayTempMin*1.2 arrayTempMax*1.2]);             
            text(  numSpatialBins*0.75, arrayTempMax*0.85, '\bfCytoplasm',  'HorizontalAlignment', 'center',  'FontSize', numPlotFontSizeTitle  );
            set(gca, 'XTick', arrayLayerBoundaries, 'XTickLabel', [0 1 2 3]);
            set(gca, 'YTick', [ -4, -3, -2, -1, 0, 1, 2, 3, 4 ], 'YTickLabel', []);
            for iTick = ceil(arrayTempMin):floor(arrayTempMax),
                if (iTick >= 2),
                    text(-1, double(iTick), ['$$\bar{\mu}$$+' num2str(iTick) '$$\bar{\sigma}$$'], 'Interpreter', 'latex', 'HorizontalAlignment', 'right');
                elseif (iTick <= -2),
                    text(-1, double(iTick), ['$$\bar{\mu}$$-' num2str(iTick) '$$\bar{\sigma}$$'], 'Interpreter', 'latex', 'HorizontalAlignment', 'right');
                elseif iTick == 0,
                    text(-1, double(iTick), '$$\bar{\mu}$$', 'Interpreter', 'latex', 'HorizontalAlignment', 'right');
                elseif iTick == 1
                    text(-1, double(iTick), ['$$\bar{\mu}$$+$$\bar{\sigma}$$'], 'Interpreter', 'latex', 'HorizontalAlignment', 'right');
                elseif iTick == -1
                    text(-1, double(iTick), '$$\bar{\mu}$$-$$-\bar{\sigma}$$', 'Interpreter', 'latex', 'HorizontalAlignment', 'right');
                end
            end
            set(gca, 'FontSize', numPlotFontSize);
            xlabel('d_n_o_r_m','FontSize',numPlotFontSizeTitle, 'Position', [14 arrayTempMin*1.35 1]);
            ylabel('I_s_,_n_o_r_m','FontSize',numPlotFontSizeTitle, 'Position', [-4.5 (arrayTempMax+arrayTempMin)/2 1]);

            if (  ( numLocalisationsForProtein == 1 ) ),
                text(-7.5,arrayTempMax*0.95,'E','FontSize',numSubFigLabelFontSize,'Color','k', 'FontWeight', 'bold');
            elseif (numLocalisationsForProtein == 2) && (length(arrayOtherImageTypes{iOutputProtein}) == 1),
                text(-7.5,arrayTempMax*0.95,'E','FontSize',numSubFigLabelFontSize,'Color','k', 'FontWeight', 'bold');
            elseif (numLocalisationsForProtein == 2) && (length(arrayOtherImageTypes{iOutputProtein}) == 2),
                text(-7.5,arrayTempMax*0.95,'F','FontSize',numSubFigLabelFontSize,'Color','k', 'FontWeight', 'bold');
            elseif numLocalisationsForProtein == 3,
                text(-7.5,arrayTempMax*0.95,'F','FontSize',numSubFigLabelFontSize,'Color','k', 'FontWeight', 'bold');
            end


            %plot the nuclear localisation (if it exists)
            if ~(arrayGroupedNodes(iOutputProtein,2) == 0),
                numNode = arrayGroupedNodes(iOutputProtein,2);


                %load the data cloud for the representative image
                stringRepImageNucDataFolder = [ stringRepImageDataFolder 'N' strFolderSep ];

                structRepImageSampleData = loadSampAnalysis(stringRepImageNucDataFolder, '.mat');

                numObjects = length(structRepImageSampleData);
                numSamplesPerObject = size(structRepImageSampleData(1).NormDist,1);
                numPixelsPerSample = size(structRepImageSampleData(1).SigInt,2);
                arraySampleValues = zeros(numObjects*numSamplesPerObject, numPixelsPerSample, 'uint8');
                arraySamplePositions = zeros(numObjects*numSamplesPerObject,1, 'double');
                for iObject = 1:numObjects,
                    for iSample = 1:numSamplesPerObject,
                        arraySamplePositions(  (iObject-1)*numSamplesPerObject + iSample  ) = structRepImageSampleData(iObject).NormDist(iSample);
                        arraySampleValues(  (iObject-1)*numSamplesPerObject + iSample, :  ) = structRepImageSampleData(iObject).SigInt(iSample,:);
                    end 
                end
                arraySamplePositions = rescaleNormDist(arraySamplePositions, numSpatialBins, arrayNormDistRatio);
                numRepImageSmpMean = mean(double(arraySampleValues(:)));
                numRepImageStDev = std(double(arraySampleValues(:)));



                handlePlotThree = subplot('Position', arrayPlot3Position);
                %overlay the lowess smoothed curves and the data cloud for the representative image
                hold on;
                plot(arraySamplePositions, (double(arraySampleValues) - numRepImageSmpMean)/numRepImageStDev, '.', 'Color', arrayDataCloudColor(arrayPatientToDisplay(iOutputProtein),:), 'MarkerSize', numPlotDataCloudMarkerSize);



                for iPatient = 1:numPatients,

                    stringNucDataFolder = [ strProcDataFolder arrayOutputProtFolders{iOutputProtein} strFolderSep 'Pat_' num2str(iPatient) strFolderSep 'N' strFolderSep ];
                    structSampleData = loadSampAnalysis(stringNucDataFolder, '.mat');
                    arrayLoessData = loadLoessCurve(stringNucDataFolder, '.mat');

                    numDataObjects = length(structSampleData);
                    numSamplesPerObject = length(structSampleData(1).SigInt);
                    arraySampleData  = zeros(numDataObjects*numSamplesPerObject,1,'uint8');
                    for iObj = 1:numDataObjects,
                        numBaseIndex = (iObj-1)*numSamplesPerObject;

                        arraySampleData((numBaseIndex+1):(numBaseIndex+iSample)) = structSampleData(iObj).SigInt(:);

                    end

                    numDataMean = mean(double(arraySampleData));
                    numDataStDev = std(double(arraySampleData));



                    plot(arrayLoessData(1,:), (arrayLoessData(2,:) - numDataMean)/numDataStDev, '.', 'Color', arrayPatientColors{iPatient}, 'MarkerSize', numPlotLowessMarkerSize);
                    arrayUpperLoessCI = [arrayLoessCIBounds{1, numNode, 1}(:,1); arrayLoessCIBounds{1, numNode, 2}(:,1); arrayLoessCIBounds{1, numNode, 3}(:,1)];
                    arrayLowerLoessCI = [arrayLoessCIBounds{1, numNode, 1}(:,2); arrayLoessCIBounds{1, numNode, 2}(:,2); arrayLoessCIBounds{1, numNode, 3}(:,2)];
                    arrayLoessCIXPosCombined = [arrayLoessCIXPos{1, numNode, 1}(:); arrayLoessCIXPos{1, numNode, 2}(:); arrayLoessCIXPos{1, numNode, 3}(:)];

                    arrayLoessX = arrayLoessData(1,:);
                    arrayLoessVals = arrayLoessData(2,:);


                    for iTissueLayer = 1:numTissueLayers,

                        numTissueLayerMinX = double(arrayDivisionIndices(iTissueLayer))-1;
                        numTissueLayerMaxX = double(arrayDivisionIndices(iTissueLayer+1))-1;

                        if iTissueLayer < numTissueLayers,
                            arrayIndexInTissueLayer = find( (arrayLoessX >= numTissueLayerMinX) & (arrayLoessX < numTissueLayerMaxX) );
                            arrayCIIndexInTissueLayer = find( (arrayLoessCIXPosCombined >= numTissueLayerMinX) & (arrayLoessCIXPosCombined < numTissueLayerMaxX) );
                        else
                            arrayIndexInTissueLayer = find( (arrayLoessX >= numTissueLayerMinX) & (arrayLoessX <= numTissueLayerMaxX) );
                            arrayCIIndexInTissueLayer = find( (arrayLoessCIXPosCombined >= numTissueLayerMinX) & (arrayLoessCIXPosCombined <= numTissueLayerMaxX) );
                        end
                        arrayXinTissueLayer = arrayLoessX(arrayIndexInTissueLayer);
                        arrayCIXinTissueLayer = arrayLoessCIXPosCombined(arrayCIIndexInTissueLayer);

                        numPlotMinX = max([arrayXinTissueLayer(1) arrayCIXinTissueLayer(1)]);
                        numPlotMaxX = min([arrayXinTissueLayer(end) arrayCIXinTissueLayer(end)]);


                        arrayXToPlot = linspace(numPlotMinX, numPlotMaxX,40);

                        arrayLoessToPlot = interp1(arrayXinTissueLayer, arrayLoessVals(arrayIndexInTissueLayer), arrayXToPlot);
                        arrayUpperLoessToPlot = interp1(arrayCIXinTissueLayer, arrayUpperLoessCI(arrayCIIndexInTissueLayer), arrayXToPlot);
                        arrayLowerLoessToPlot = interp1(arrayCIXinTissueLayer, arrayLowerLoessCI(arrayCIIndexInTissueLayer), arrayXToPlot);


                        plot(arrayXToPlot, (arrayLoessToPlot - numDataMean + arrayUpperLoessToPlot)/numDataStDev, '--', 'Color', arrayPatientColors{iPatient},'LineWidth',(numPlotLowessMarkerSize/8));
                        plot(arrayXToPlot, (arrayLoessToPlot - numDataMean + arrayLowerLoessToPlot)/numDataStDev, '--', 'Color', arrayPatientColors{iPatient},'LineWidth',(numPlotLowessMarkerSize/8));

                    end
                end

                hold off;


                %extract the max and minimum values for adjusting axes
                arrayTempMax = (double(max(arraySampleValues(:))) - numRepImageSmpMean)/numRepImageStDev;
                arrayTempMin = (double(min(arraySampleValues(:))) - numRepImageSmpMean)/numRepImageStDev;     
                if arrayTempMax > 3 || isnan(arrayTempMax),
                    arrayTempMax = 3;
                end
                if arrayTempMin < -3 || isnan(arrayTempMin),
                    arrayTempMin = -3;
                end

                %set the axes
                axis([0 numSpatialBins, arrayTempMin*1.2 arrayTempMax*1.2]);
                text(  numSpatialBins*0.75, arrayTempMax*0.85, '\bfNucleus',  'HorizontalAlignment', 'center',  'FontSize', numPlotFontSizeTitle  );
                set(gca, 'XTick', arrayLayerBoundaries, 'XTickLabel', [0 1 2 3]);
                set(gca, 'YTick', [ -4, -3, -2, -1, 0, 1, 2, 3, 4 ], 'YTickLabel', []);
                for iTick = ceil(arrayTempMin):floor(arrayTempMax),
                    if (iTick >= 2),
                        text(-1, double(iTick), ['$$\bar{\mu}$$+' num2str(iTick) '$$\bar{\sigma}$$'], 'Interpreter', 'latex', 'HorizontalAlignment', 'right');
                    elseif (iTick <= -2),
                        text(-1, double(iTick), ['$$\bar{\mu}$$-' num2str(iTick) '$$\bar{\sigma}$$'], 'Interpreter', 'latex', 'HorizontalAlignment', 'right');
                    elseif iTick == 0,
                        text(-1, double(iTick), '$$\bar{\mu}$$', 'Interpreter', 'latex', 'HorizontalAlignment', 'right');
                    elseif iTick == 1
                        text(-1, double(iTick), ['$$\bar{\mu}$$+$$\bar{\sigma}$$'], 'Interpreter', 'latex', 'HorizontalAlignment', 'right');
                    elseif iTick == -1
                        text(-1, double(iTick), '$$\bar{\mu}$$-$$-\bar{\sigma}$$', 'Interpreter', 'latex', 'HorizontalAlignment', 'right');
                    end
                end
                set(gca, 'FontSize', numPlotFontSize);
                xlabel('d_n_o_r_m','FontSize',numPlotFontSizeTitle, 'Position', [14 arrayTempMin*1.35 1]);
                ylabel('I_s_,_n_o_r_m','FontSize',numPlotFontSizeTitle, 'Position', [-4.5 (arrayTempMax+arrayTempMin)/2 1]);


                text(-7.5,arrayTempMax*0.95,'G','FontSize',numSubFigLabelFontSize,'Color','k', 'FontWeight', 'bold');

            end

        end

    else
        %the cytoplasmic localisation is always sampled, plot it first
        if ~(arrayGroupedNodes(iOutputProtein,1) == 0),
            numNode = arrayGroupedNodes(iOutputProtein,1);

            %load the data cloud for the representative image
            stringRepImageCytoDataFolder = [ stringRepImageDataFolder 'C' strFolderSep ];

            structRepImageSampleData = loadSampAnalysis(stringRepImageCytoDataFolder, '.mat');

            numObjects = length(structRepImageSampleData);
            numSamplesPerObject = size(structRepImageSampleData(1).NormDist,1);
            numPixelsPerSample = size(structRepImageSampleData(1).SigInt,2);
            arraySampleValues = zeros(numObjects*numSamplesPerObject, numPixelsPerSample, 'uint8');
            arraySamplePositions = zeros(numObjects*numSamplesPerObject,1, 'double');
            for iObject = 1:numObjects,
                for iSample = 1:numSamplesPerObject,
                    arraySamplePositions(  (iObject-1)*numSamplesPerObject + iSample  ) = structRepImageSampleData(iObject).NormDist(iSample);
                    arraySampleValues(  (iObject-1)*numSamplesPerObject + iSample, :  ) = structRepImageSampleData(iObject).SigInt(iSample,:);
                end 
            end
            arraySamplePositions = rescaleNormDist(arraySamplePositions, numSpatialBins, arrayNormDistRatio);
            numRepImageSmpMean = mean(double(arraySampleValues(:)));
            numRepImageStDev = std(double(arraySampleValues(:)));
            
            handlePlotOne = subplot('Position', arrayPlot1Position);

            %overlay the lowess smoothed curves and the data cloud for the representative image
            hold on;
            plot(arraySamplePositions, (double(arraySampleValues) - numRepImageSmpMean)/numRepImageStDev, '.', 'Color', arrayDataCloudColor(arrayPatientToDisplay(iOutputProtein),:), 'MarkerSize', numPlotDataCloudMarkerSize);


            for iPatient = 1:numPatients,

                stringCytoDataFolder = [ strProcDataFolder arrayOutputProtFolders{iOutputProtein} strFolderSep 'Pat_' num2str(iPatient) strFolderSep 'C' strFolderSep ];
                structSampleData = loadSampAnalysis(stringCytoDataFolder, '.mat');
                arrayLoessData = loadLoessCurve(stringCytoDataFolder, '.mat');

                numDataObjects = length(structSampleData);
                numSamplesPerObject = length(structSampleData(1).SigInt);
                arraySampleData  = zeros(numDataObjects*numSamplesPerObject,1,'uint8');
                for iObj = 1:numDataObjects,
                    numBaseIndex = (iObj-1)*numSamplesPerObject;

                    arraySampleData((numBaseIndex+1):(numBaseIndex+iSample)) = structSampleData(iObj).SigInt(:);

                end

                numDataMean = mean(double(arraySampleData));
                numDataStDev = std(double(arraySampleData));

                plot(arrayLoessData(1,:), (arrayLoessData(2,:) - numDataMean)/numDataStDev, '.', 'Color', arrayPatientColors{iPatient}, 'MarkerSize', numPlotLowessMarkerSize);
                arrayUpperLoessCI = [arrayLoessCIBounds{1, numNode, 1}(:,1); arrayLoessCIBounds{1, numNode, 2}(:,1); arrayLoessCIBounds{1, numNode, 3}(:,1)];
                arrayLowerLoessCI = [arrayLoessCIBounds{1, numNode, 1}(:,2); arrayLoessCIBounds{1, numNode, 2}(:,2); arrayLoessCIBounds{1, numNode, 3}(:,2)];

                arrayLoessCIXPosCombined = [arrayLoessCIXPos{1, numNode, 1}(:); arrayLoessCIXPos{1, numNode, 2}(:); arrayLoessCIXPos{1, numNode, 3}(:)];

                arrayLoessX = arrayLoessData(1,:);
                arrayLoessVals = arrayLoessData(2,:);


                for iTissueLayer = 1:numTissueLayers,

                    numTissueLayerMinX = double(arrayDivisionIndices(iTissueLayer))-1;
                    numTissueLayerMaxX = double(arrayDivisionIndices(iTissueLayer+1))-1;

                    if iTissueLayer < numTissueLayers,
                        arrayIndexInTissueLayer = find( (arrayLoessX >= numTissueLayerMinX) & (arrayLoessX < numTissueLayerMaxX) );
                        arrayCIIndexInTissueLayer = find( (arrayLoessCIXPosCombined >= numTissueLayerMinX) & (arrayLoessCIXPosCombined < numTissueLayerMaxX) );
                    else
                        arrayIndexInTissueLayer = find( (arrayLoessX >= numTissueLayerMinX) & (arrayLoessX <= numTissueLayerMaxX) );
                        arrayCIIndexInTissueLayer = find( (arrayLoessCIXPosCombined >= numTissueLayerMinX) & (arrayLoessCIXPosCombined <= numTissueLayerMaxX) );
                    end
                    arrayXinTissueLayer = arrayLoessX(arrayIndexInTissueLayer);
                    arrayCIXinTissueLayer = arrayLoessCIXPosCombined(arrayCIIndexInTissueLayer);

                    numPlotMinX = max([arrayXinTissueLayer(1) arrayCIXinTissueLayer(1)]);
                    numPlotMaxX = min([arrayXinTissueLayer(end) arrayCIXinTissueLayer(end)]);


                    arrayXToPlot = linspace(numPlotMinX, numPlotMaxX,40);
                    arrayLoessToPlot = interp1(arrayXinTissueLayer, arrayLoessVals(arrayIndexInTissueLayer), arrayXToPlot);
                    arrayUpperLoessToPlot = interp1(arrayCIXinTissueLayer, arrayUpperLoessCI(arrayCIIndexInTissueLayer), arrayXToPlot);
                    arrayLowerLoessToPlot = interp1(arrayCIXinTissueLayer, arrayLowerLoessCI(arrayCIIndexInTissueLayer), arrayXToPlot);


                    plot(arrayXToPlot, (arrayLoessToPlot - numDataMean + arrayUpperLoessToPlot)/numDataStDev, '--', 'Color', arrayPatientColors{iPatient},'LineWidth',(numPlotLowessMarkerSize/8));
                    plot(arrayXToPlot, (arrayLoessToPlot - numDataMean + arrayLowerLoessToPlot)/numDataStDev, '--', 'Color', arrayPatientColors{iPatient},'LineWidth',(numPlotLowessMarkerSize/8));

                end
            end

            hold off;

            %extract the max and minimum values for adjusting axes
            arrayTempMax = (double(max(arraySampleValues(:))) - numRepImageSmpMean)/numRepImageStDev;
            arrayTempMin = (double(min(arraySampleValues(:))) - numRepImageSmpMean)/numRepImageStDev;     
            if arrayTempMax > 3 || isnan(arrayTempMax),
                arrayTempMax = 3;
            end
            if arrayTempMin < -3 || isnan(arrayTempMin),
                arrayTempMin = -3;
            end

            %set the axes
            axis([0 numSpatialBins, arrayTempMin*1.2 arrayTempMax*1.2]);
            text(  numSpatialBins*0.75, arrayTempMax*0.85, '\bfCytoplasm',  'HorizontalAlignment', 'center',  'FontSize', numPlotFontSizeTitle  );
            set(gca, 'XTick', arrayLayerBoundaries, 'XTickLabel', [0 1 2 3]);
            set(gca, 'YTick', [ -4, -3, -2, -1, 0, 1, 2, 3, 4 ], 'YTickLabel', []);
            for iTick = ceil(arrayTempMin):floor(arrayTempMax),
                if (iTick >= 2),
                    text(-1, double(iTick), ['$$\bar{\mu}$$+' num2str(iTick) '$$\bar{\sigma}$$'], 'Interpreter', 'latex', 'HorizontalAlignment', 'right');
                elseif (iTick <= -2),
                    text(-1, double(iTick), ['$$\bar{\mu}$$-' num2str(iTick) '$$\bar{\sigma}$$'], 'Interpreter', 'latex', 'HorizontalAlignment', 'right');
                elseif iTick == 0,
                    text(-1, double(iTick), '$$\bar{\mu}$$', 'Interpreter', 'latex', 'HorizontalAlignment', 'right');
                elseif iTick == 1
                    text(-1, double(iTick), ['$$\bar{\mu}$$+$$\bar{\sigma}$$'], 'Interpreter', 'latex', 'HorizontalAlignment', 'right');
                elseif iTick == -1
                    text(-1, double(iTick), '$$\bar{\mu}$$-$$-\bar{\sigma}$$', 'Interpreter', 'latex', 'HorizontalAlignment', 'right');
                end
            end
            set(gca, 'FontSize', numPlotFontSize);
            xlabel('d_n_o_r_m','FontSize',numPlotFontSizeTitle, 'Position', [14 arrayTempMin*1.35 1]);
            ylabel('I_s_,_n_o_r_m','FontSize',numPlotFontSizeTitle, 'Position', [-4.5 (arrayTempMax+arrayTempMin)/2 1]);


            if (  ( numLocalisationsForProtein == 1 ) ),
                text(-7.5,arrayTempMax*0.95,'D','FontSize',numSubFigLabelFontSize,'Color','k', 'FontWeight', 'bold');
            elseif (numLocalisationsForProtein == 2) && (length(arrayOtherImageTypes{iOutputProtein}) == 1),
                text(-7.5,arrayTempMax*0.95,'D','FontSize',numSubFigLabelFontSize,'Color','k', 'FontWeight', 'bold');
            elseif (numLocalisationsForProtein == 2) && (length(arrayOtherImageTypes{iOutputProtein}) == 2),
                text(-7.5,arrayTempMax*0.95,'E','FontSize',numSubFigLabelFontSize,'Color','k', 'FontWeight', 'bold');
            elseif numLocalisationsForProtein == 3,
                text(-7.5,arrayTempMax*0.95,'E','FontSize',numSubFigLabelFontSize,'Color','k', 'FontWeight', 'bold');
            end


            %plot the nuclear localisation (if it exists)
            if ~(arrayGroupedNodes(iOutputProtein,2) == 0),
                numNode = arrayGroupedNodes(iOutputProtein,2);

                %load the data cloud for the representative image
                stringRepImageNucDataFolder = [ stringRepImageDataFolder 'N' strFolderSep ];

                structRepImageSampleData = loadSampAnalysis(stringRepImageNucDataFolder, '.mat');

                numObjects = length(structRepImageSampleData);
                numSamplesPerObject = size(structRepImageSampleData(1).NormDist,1);
                numPixelsPerSample = size(structRepImageSampleData(1).SigInt,2);
                arraySampleValues = zeros(numObjects*numSamplesPerObject, numPixelsPerSample, 'uint8');
                arraySamplePositions = zeros(numObjects*numSamplesPerObject,1, 'double');
                for iObject = 1:numObjects,
                    for iSample = 1:numSamplesPerObject,
                        arraySamplePositions(  (iObject-1)*numSamplesPerObject + iSample  ) = structRepImageSampleData(iObject).NormDist(iSample);
                        arraySampleValues(  (iObject-1)*numSamplesPerObject + iSample, :  ) = structRepImageSampleData(iObject).SigInt(iSample,:);
                    end 
                end
                arraySamplePositions = rescaleNormDist(arraySamplePositions, numSpatialBins, arrayNormDistRatio);
                numRepImageSmpMean = mean(double(arraySampleValues(:)));
                numRepImageStDev = std(double(arraySampleValues(:)));


                handlePlotTwo = subplot('Position', arrayPlot2Position);


                %overlay the lowess smoothed curves and the data cloud for the representative image
                hold on;
                plot(arraySamplePositions, (double(arraySampleValues) - numRepImageSmpMean)/numRepImageStDev, '.', 'Color', arrayDataCloudColor(arrayPatientToDisplay(iOutputProtein),:), 'MarkerSize', numPlotDataCloudMarkerSize);


                for iPatient = 1:numPatients,

                    stringNucDataFolder = [ strProcDataFolder arrayOutputProtFolders{iOutputProtein} strFolderSep 'Pat_' num2str(iPatient) strFolderSep 'N' strFolderSep ];
                    structSampleData = loadSampAnalysis(stringNucDataFolder, '.mat');
                    arrayLoessData = loadLoessCurve(stringNucDataFolder, '.mat');

                    numDataObjects = length(structSampleData);
                    numSamplesPerObject = length(structSampleData(1).SigInt);
                    arraySampleData  = zeros(numDataObjects*numSamplesPerObject,1,'uint8');
                    for iObj = 1:numDataObjects,
                        numBaseIndex = (iObj-1)*numSamplesPerObject;

                        arraySampleData((numBaseIndex+1):(numBaseIndex+iSample)) = structSampleData(iObj).SigInt(:);

                    end

                    numDataMean = mean(double(arraySampleData));
                    numDataStDev = std(double(arraySampleData));



                    plot(arrayLoessData(1,:), (arrayLoessData(2,:) - numDataMean)/numDataStDev, '.', 'Color', arrayPatientColors{iPatient}, 'MarkerSize', numPlotLowessMarkerSize);
                    arrayUpperLoessCI = [arrayLoessCIBounds{1, numNode, 1}(:,1); arrayLoessCIBounds{1, numNode, 2}(:,1); arrayLoessCIBounds{1, numNode, 3}(:,1)];
                    arrayLowerLoessCI = [arrayLoessCIBounds{1, numNode, 1}(:,2); arrayLoessCIBounds{1, numNode, 2}(:,2); arrayLoessCIBounds{1, numNode, 3}(:,2)];

                    arrayLoessCIXPosCombined = [arrayLoessCIXPos{1, numNode, 1}(:); arrayLoessCIXPos{1, numNode, 2}(:); arrayLoessCIXPos{1, numNode, 3}(:)];

                    arrayLoessX = arrayLoessData(1,:);
                    arrayLoessVals = arrayLoessData(2,:);


                    for iTissueLayer = 1:numTissueLayers,

                        numTissueLayerMinX = double(arrayDivisionIndices(iTissueLayer))-1;
                        numTissueLayerMaxX = double(arrayDivisionIndices(iTissueLayer+1))-1;

                        if iTissueLayer < numTissueLayers,
                            arrayIndexInTissueLayer = find( (arrayLoessX >= numTissueLayerMinX) & (arrayLoessX < numTissueLayerMaxX) );
                            arrayCIIndexInTissueLayer = find( (arrayLoessCIXPosCombined >= numTissueLayerMinX) & (arrayLoessCIXPosCombined < numTissueLayerMaxX) );
                        else
                            arrayIndexInTissueLayer = find( (arrayLoessX >= numTissueLayerMinX) & (arrayLoessX <= numTissueLayerMaxX) );
                            arrayCIIndexInTissueLayer = find( (arrayLoessCIXPosCombined >= numTissueLayerMinX) & (arrayLoessCIXPosCombined <= numTissueLayerMaxX) );
                        end
                        arrayXinTissueLayer = arrayLoessX(arrayIndexInTissueLayer);
                        arrayCIXinTissueLayer = arrayLoessCIXPosCombined(arrayCIIndexInTissueLayer);

                        numPlotMinX = max([arrayXinTissueLayer(1) arrayCIXinTissueLayer(1)]);
                        numPlotMaxX = min([arrayXinTissueLayer(end) arrayCIXinTissueLayer(end)]);


                        arrayXToPlot = linspace(numPlotMinX, numPlotMaxX,40);

                        arrayLoessToPlot = interp1(arrayXinTissueLayer, arrayLoessVals(arrayIndexInTissueLayer), arrayXToPlot);
                        arrayUpperLoessToPlot = interp1(arrayCIXinTissueLayer, arrayUpperLoessCI(arrayCIIndexInTissueLayer), arrayXToPlot);
                        arrayLowerLoessToPlot = interp1(arrayCIXinTissueLayer, arrayLowerLoessCI(arrayCIIndexInTissueLayer), arrayXToPlot);


                        plot(arrayXToPlot, (arrayLoessToPlot - numDataMean + arrayUpperLoessToPlot)/numDataStDev, '--', 'Color', arrayPatientColors{iPatient},'LineWidth',(numPlotLowessMarkerSize/8));
                        plot(arrayXToPlot, (arrayLoessToPlot - numDataMean + arrayLowerLoessToPlot)/numDataStDev, '--', 'Color', arrayPatientColors{iPatient},'LineWidth',(numPlotLowessMarkerSize/8));

                    end
                end

                hold off;


                %extract the max and minimum values for adjusting axes
                arrayTempMax = (double(max(arraySampleValues(:))) - numRepImageSmpMean)/numRepImageStDev;
                arrayTempMin = (double(min(arraySampleValues(:))) - numRepImageSmpMean)/numRepImageStDev;     
                if arrayTempMax > 3 || isnan(arrayTempMax),
                    arrayTempMax = 3;
                end
                if arrayTempMin < -3 || isnan(arrayTempMin),
                    arrayTempMin = -3;
                end

                %set the axes
                axis([0 numSpatialBins, arrayTempMin*1.2 arrayTempMax*1.2]);
                text(  numSpatialBins*0.75, arrayTempMax*0.85, '\bfNucleus',  'HorizontalAlignment', 'center',  'FontSize', numPlotFontSizeTitle  );
                set(gca, 'XTick', arrayLayerBoundaries, 'XTickLabel', [0 1 2 3]);
                set(gca, 'YTick', [ -4, -3, -2, -1, 0, 1, 2, 3, 4 ], 'YTickLabel', []);
                for iTick = ceil(arrayTempMin):floor(arrayTempMax),
                    if (iTick >= 2),
                        text(-1, double(iTick), ['$$\bar{\mu}$$+' num2str(iTick) '$$\bar{\sigma}$$'], 'Interpreter', 'latex', 'HorizontalAlignment', 'right');
                    elseif (iTick <= -2),
                        text(-1, double(iTick), ['$$\bar{\mu}$$-' num2str(iTick) '$$\bar{\sigma}$$'], 'Interpreter', 'latex', 'HorizontalAlignment', 'right');
                    elseif iTick == 0,
                        text(-1, double(iTick), '$$\bar{\mu}$$', 'Interpreter', 'latex', 'HorizontalAlignment', 'right');
                    elseif iTick == 1
                        text(-1, double(iTick), ['$$\bar{\mu}$$+$$\bar{\sigma}$$'], 'Interpreter', 'latex', 'HorizontalAlignment', 'right');
                    elseif iTick == -1
                        text(-1, double(iTick), '$$\bar{\mu}$$-$$-\bar{\sigma}$$', 'Interpreter', 'latex', 'HorizontalAlignment', 'right');
                    end
                end
                set(gca, 'FontSize', numPlotFontSize);
                xlabel('d_n_o_r_m','FontSize',numPlotFontSizeTitle, 'Position', [14 arrayTempMin*1.35 1]);
                ylabel('I_s_,_n_o_r_m','FontSize',numPlotFontSizeTitle, 'Position', [-4.5 (arrayTempMax+arrayTempMin)/2 1]);

                if (  ( numLocalisationsForProtein == 1 ) ),
                    text(-7.5,arrayTempMax*0.95,'E','FontSize',numSubFigLabelFontSize,'Color','k', 'FontWeight', 'bold');
                elseif (numLocalisationsForProtein == 2) && (length(arrayOtherImageTypes{iOutputProtein}) == 1),
                    text(-7.5,arrayTempMax*0.95,'E','FontSize',numSubFigLabelFontSize,'Color','k', 'FontWeight', 'bold');
                elseif (numLocalisationsForProtein == 2) && (length(arrayOtherImageTypes{iOutputProtein}) == 2),
                    text(-7.5,arrayTempMax*0.95,'F','FontSize',numSubFigLabelFontSize,'Color','k', 'FontWeight', 'bold');
                elseif numLocalisationsForProtein == 3,
                    text(-7.5,arrayTempMax*0.95,'F','FontSize',numSubFigLabelFontSize,'Color','k', 'FontWeight', 'bold');
                end


            end

        end

    end
%     end
    
    if (iOutputProtein >= 7) && (iOutputProtein <= 10),
    %MEK, pMEK, ERK or pERK; plot the cyto:nuc signal ratio in the
    % final panel
        handleCytoNucRatio = subplot('Position', arrayPlot3Position);
        
        stringProteinDataFolder = [ strProcDataFolder arrayOutputProtFolders{iOutputProtein} strFolderSep ];
       
        hold on;
        
        for iPatient = 1:numPatients,
            
            stringCytoDataFolder = [ strProcDataFolder arrayOutputProtFolders{iOutputProtein} strFolderSep 'Pat_' num2str(iPatient) strFolderSep 'C' strFolderSep ];
            stringNucDataFolder = [ strProcDataFolder arrayOutputProtFolders{iOutputProtein} strFolderSep 'Pat_' num2str(iPatient) strFolderSep 'N' strFolderSep ];
            arrayLoessCytoData = loadLoessCurve(stringCytoDataFolder, '.mat');
            arrayLoessNucData = loadLoessCurve(stringNucDataFolder, '.mat');
            
            
            for iTissueLayer = 1:numTissueLayers,

                numTissueLayerMinX = double(arrayDivisionIndices(iTissueLayer))-1;
                numTissueLayerMaxX = double(arrayDivisionIndices(iTissueLayer+1))-1;

                if iTissueLayer < numTissueLayers,
                    arrayCytoIndexInTissueLayer = find( (arrayLoessCytoData(1,:) >= numTissueLayerMinX) & (arrayLoessCytoData(1,:) < numTissueLayerMaxX) );
                    arrayNucIndexInTissueLayer = find( (arrayLoessNucData(1,:) >= numTissueLayerMinX) & (arrayLoessNucData(1,:) < numTissueLayerMaxX) );
                else
                    arrayCytoIndexInTissueLayer = find( (arrayLoessCytoData(1,:) >= numTissueLayerMinX) & (arrayLoessCytoData(1,:) <= numTissueLayerMaxX) );
                    arrayNucIndexInTissueLayer = find( (arrayLoessNucData(1,:) >= numTissueLayerMinX) & (arrayLoessNucData(1,:) <= numTissueLayerMaxX) );
                end
            
                numStartX = max([ arrayLoessCytoData(1,arrayCytoIndexInTissueLayer(1)) arrayLoessNucData(1,arrayNucIndexInTissueLayer(1)) ]);
                numEndX = min([ arrayLoessCytoData(1,arrayCytoIndexInTissueLayer(end)) arrayLoessNucData(1,arrayNucIndexInTissueLayer(end)) ]);

                arrayInterpXPos = linspace(numStartX, numEndX, numSpatialBins*4);

                arrayInterpNucData = interp1(arrayLoessNucData(1,:), arrayLoessNucData(2,:), arrayInterpXPos);
                arrayInterpCytoData = interp1(arrayLoessCytoData(1,:), arrayLoessCytoData(2,:), arrayInterpXPos);

                plot(arrayInterpXPos, arrayInterpCytoData./arrayInterpNucData, '-', 'Color', arrayPatientColors{iPatient,:}, 'LineWidth', 3);
                
            end
        end
              
        hold off;
        
        arrayYLim = get(gca, 'YLim');
        arrayTempMin = arrayYLim(1);
        arrayTempMax = arrayYLim(2);
        
        if arrayTempMin > 0.5,
            arrayTempMin = 0.5;
        end
        if arrayTempMax < 2,
            arrayTempMax = 2;
        end
        set(gca, 'YLim', [arrayTempMin arrayTempMax]);
        
        %set the axes
        set(gca, 'XLim', [0 numSpatialBins]);
        set(gca, 'XTick', arrayLayerBoundaries, 'XTickLabel', [0 1 2 3]);
        set(gca, 'FontSize', numPlotFontSize);
        xlabel('d_n_o_r_m','FontSize',numPlotFontSizeTitle, 'Position', [14 arrayTempMin-(arrayTempMax - arrayTempMin)*0.1 1]);
        ylabel('(I_s_,_c_y_t_o)/(I_s_,_n_u_c)','FontSize',numPlotFontSizeTitle, 'Position', [-4.5 (arrayTempMax+arrayTempMin)/2 1]);
        text(-7.5,arrayTempMax*0.95,'F','FontSize',numSubFigLabelFontSize,'Color','k', 'FontWeight', 'bold');
        
    end
                     
    print(figOut, '-r300', '-dpng', [stringOutputDataFolder arrayOutputProtFolders{iOutputProtein} '.png']);
    close(figOut);
        
end