%% examine_epidermal_thickness.m
% This MATLAB script reproduces Fig. AF1.3 from the GigaScience Data Note:
%   Cursons et al. (2015). Spatially-transformed fluorescence image data 
%    for ERK-MAPK and selected proteins within human epidermis.
%    GigaScience. Accepted Nov 2015.
%    doi: not-yet-known
%
% It has been written to operate together with immunofluorescence image
%  data of human epidermis hosted on the GigaScience Database:
%   Cursons, J; Angel, C, E; Hurley, D, G; Print, C, G; Dunbar, P; 
%    Jacobs, M, D; Crampin, E, J (2015): Supporting data for 
%    "Spatially-transformed fluorescence image data for ERK-MAPK and
%    selected proteins within human epidermis". GigaScience Database. 
%    http://dx.doi.org/10.5524/100168
%
% A number of functions are used by this script, some of which have
%  dependencies upon MATLAB Toolboxes:
%   - Image Processing Toolbox: this script calls on a number of functions
%                                   which use the imread function
%   - Neural Network Toolbox: this script directly calls the dist function
%   - Statistics and Machine Learning Toolbox: this script directly calls
%                                   the ksdensity function
%
% This script was created by Joe Cursons at the University of Melbourne
%   Systems Biology Laboratory:
%       joseph.cursons@unimelb.edu.au
%
% Last Updated: 30/11/15
%
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Input Parameters
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  

%move through all patients contained within this data set
arrayPatients = 1:3;
numPatients = length(arrayPatients);

%move through all IF targets examined within these data
arrayTargetPaths = { ['CALM'];
                     ['ERK'];
                     ['ERK_ph'];
                     ['FOSC'];
                     ['FRA2'];
                     ['ITGB1'];
                     ['ITGB4'];
                     ['JUNB'];
                     ['JUNC'];
                     ['K10'];
                     ['K14'];
                     ['MEK'];
                     ['MEK_ph'];
                     ['RAF'];
                     ['RAF_ph'];
                     ['SFN'] };
numTargets = length(arrayTargetPaths);


%measure the thickness of each segmented tissue layer (basal, spinous &
% granular, transitional) and the complete epidermis
numDistToMeasure = 4;

%specify strings for the distances being measured 
arrayDistStrings = { {['Basal'];['layer']};
                     {['Spinous and'];['granular layers']};
                     {['Transitional'];['layer']};
                     {['Total'];['epidermis']} };

%specify strings for labelling the sub-figure panels
arraySubFigStrings = { ['A'];
                       ['B'];
                       ['C'];
                       ['D'] };

%define plotting colors for the patients:
arrayPatientColors = { [ 1 0 0 ];   %   Pat 1 - Red
                       [ 0 1 0 ];   %   Pat 2 - Green
                       [ 0 0 1 ] }; %   Pat 3 - Blue

%define the image pixel resolutions (um per pixel)
arrayPixRes = [ 0.135444, 0.132890, 0.128008;       %CALM
                0.147479, 0.159570, 0.149920;       %ERK
                0.145719, 0.232515, 0.151396;       %ERK_ph
                0.149750, 0.156051, 0.110638;       %FOSC
                0.144868, 0.163601, 0.118187;       %FRA2
                0.148160, 0.157981, 0.144300;       %ITGB1
                0.142824, 0.154177, 0.170356;       %ITGB4
                0.140951, 0.141575, 0.129257;       %JUNB
                0.138680, 0.154121, 0.129257;       %JUNC
                0.115122, 0.132266, 0.149012;       %K10
                0.125964, 0.137829, 0.126135;       %K14
                0.139532, 0.134139, 0.119153;       %MEK
                0.142200, 0.172343, 0.120401;       %MEK_ph
                0.168028, 0.173308, 0.154915;       %RAF
                0.105302, 0.137147, 0.162636;       %RAF_ph
                0.152020, 0.144811, 0.160251 ];     %SFN
                      

%specify output plotting settings
arrayOutputPlotPaperSize = [ 20, 20 ];
arrayOutputPlotPaperPosition = [ 0, 0, 20, 20 ];
arrayOutputPlotScreenPos = [ 50 50 1050 1050 ];

 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Path Manipulations
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  

%depending upon the OS, use a different folder separater (forward vs back
% slash)
if ispc,
    strFoldSep = '\';
elseif isunix,
    strFoldSep = '/';
else
    disp(['warning: cannot determine the operating system, defaulting to ' ...
            'forward slash for the folder path separator' ] );
    strFoldSep = '/';
end
 
%determine the current directory
strCurrDir = cd;

%the functions used by this script should be located in the /functions/
% folder, so add this to the MATLAB file path
if ~isdeployed
    addpath(genpath(strCurrDir));
else
    addpath(genpath([ctfroot '/code']))
end


%manipulate the file path to determine the appropriate folders
arrayCurrDirFoldSepPos = strfind(strCurrDir, strFoldSep);

if arrayCurrDirFoldSepPos(end) == length(strCurrDir),
    %there is a backslash at the end
    strBaseDir = strCurrDir(1:(arrayCurrDirFoldSepPos(end-1)));
else
    strBaseDir = strCurrDir(1:(arrayCurrDirFoldSepPos(end)));
end

%and just output to the root directory
stringOutputDir = strBaseDir;


%and the relative path of the data
arrayCurrDirSeps = strfind(strCurrDir, strFoldSep);
numDirSeps = length(arrayCurrDirSeps);
stringMountedFolder = strCurrDir(1:arrayCurrDirSeps(numDirSeps));


 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Load the layer boundary data and calculate minimum distances
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
 
%create an output array for each measured thickness
arrayLayerPixThickness = cell(numPatients, numTargets, numDistToMeasure);

%create an output array for the number of z-slices sampled for each image
% data stack
arrayZSlices = cell(numPatients,numTargets);

%move through all specified targets
disp([char(9) 'calculating the thickness of each segmented epidermal tissue layer and the']);
disp([char(9) '  full epidermis, this may take some time..']);
for iTarget = 1:numTargets,    
    
    numPercComplete = 100.0*double(iTarget-1)/double(numTargets);
    disp([char(9) char(9) 'calculating for: ' arrayTargetPaths{iTarget} '; approx. ' num2str(numPercComplete, '%03.1f') '% complete']);

    %determine the relative image data folder path
    strImgDataFolder = [ stringMountedFolder 'image' strFoldSep arrayTargetPaths{iTarget} strFoldSep ];

    %determine the relative processed data folder path
    strProcImgDataPath = [ stringMountedFolder 'processed' strFoldSep arrayTargetPaths{iTarget} strFoldSep ];
    
    %and move through all patients
    for iPatient = 1:numPatients,

        %extract the patient number
        numPatient = arrayPatients(iPatient);
        %and then use this to specify the data folder
        stringSegImgPath = [ strProcImgDataPath 'Pat_' num2str(iPatient) strFoldSep ];
        %use the cytoplasmic sample location dat (exists for all targets)
        stringSmpLocPath = [ strProcImgDataPath 'Pat_' num2str(iPatient) strFoldSep 'C' strFoldSep ];
        %specify the original image data folder for this patient
        strPatImgDataFolder = [ strImgDataFolder 'Pat_' num2str(iPatient) strFoldSep 'image_data' strFoldSep ];

        %check that the specified sub-cellular localisation exists
        % for processing
        if (exist(stringSegImgPath, 'dir') == 7) && (exist(stringSmpLocPath, 'dir') == 7),

            %load the sample location data, as it is the cleanest
            % way to determine sampled z-positions for loading the 
            % sparse cell arrays
            structSampleLocData = loadSampLocs(stringSmpLocPath, '.mat');

            %examine the sample data locations to determine the z-positions at 
            % which samples have been collected
            arrayAllSliceZPositions = zeros(length(structSampleLocData),1,'uint8');
            for iSampleGroup = 1:length(structSampleLocData),
                arrayAllSliceZPositions(iSampleGroup) = structSampleLocData(iSampleGroup).ZPosition;
            end
            %make this list unique
            arrayZSlicesSampled = unique(arrayAllSliceZPositions);
            
            %save to the output array
            arrayZSlices{iPatient,iTarget} = arrayZSlicesSampled;
            
            %load the full set of tissue segmentation data
            [ CellImageStack, CellImageLayers, CellImageBasalLamina, CellImageBoundaryOne, CellImageBoundaryTwo, CellImageOuterBoundary ] = ...
                loadImageStackAsSparse3DCellArrays( arrayZSlicesSampled, strPatImgDataFolder, stringSegImgPath);

            %specify the size of the final output array for the distances
            numZSlices = length(arrayZSlices{iPatient,iTarget});
            for iDist = 1:numDistToMeasure,
                arrayLayerPixThickness{iPatient, iTarget, iDist} = cell(numZSlices);
            end
            
            %move through each sampled/segmented z-slice
            for iZPos = 1:numZSlices,
                
                %extract the index for this z-slice within the data
                numZPos = arrayZSlices{iPatient,iTarget}(iZPos);
                
                %calculate each distance
                for iDist = 1:numDistToMeasure,
                    
                    %load the masked boundary images corresponding to the
                    % distance desired
                    if iDist == 1,
                        imageLowerBound = CellImageBasalLamina{numZPos};
                        imageUpperBound = CellImageBoundaryOne{numZPos};
                    elseif iDist == 2,
                        imageLowerBound = CellImageBoundaryOne{numZPos};
                        imageUpperBound = CellImageBoundaryTwo{numZPos};
                    elseif iDist == 3,
                        imageLowerBound = CellImageBoundaryTwo{numZPos};
                        imageUpperBound = CellImageOuterBoundary{numZPos};
                    elseif iDist == 4,
                        imageLowerBound = CellImageBasalLamina{numZPos};
                        imageUpperBound = CellImageOuterBoundary{numZPos};
                    end
                    
                    %extract the pixel co-ordinates
                    [arrayLowerY, arrayLowerX] = find(imageLowerBound);
                    [arrayUpperY, arrayUpperX] = find(imageUpperBound);
                    numLowerPix = length(arrayLowerY);
                    numUpperPix = length(arrayUpperY);

                    %create the output vector
                    arrayLayerPixThickness{iPatient, iTarget, iDist}{iZPos} = zeros(numLowerPix,1,'double');

                    %combine the pixel lists and use the dist function
                    arrayCombinedPix = [ arrayLowerX', arrayUpperX'; arrayLowerY', arrayUpperY' ];
                    arrayFullDistMatrix = dist(arrayCombinedPix);

                    %extract the section of the distance matrix which is of interest
                    arrayLowerToUpper = arrayFullDistMatrix(1:numLowerPix,numLowerPix+1:end);

                    %and determine the minimum values
                    for iPix = 1:numLowerPix,
                        arrayLayerPixThickness{iPatient, iTarget, iDist}{iZPos}(iPix) = min(arrayLowerToUpper(iPix,:));
                    end

                end
            end
            
        end
              
    end
end

 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Perform post-processing to combine the distances across all targets and 
%   rescale to absolute distance
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
arrayCombAllObs = cell(numPatients,numDistToMeasure);
for iPatient = 1:numPatients,
    for iDist = 1:numDistToMeasure,
        
        %create the combined/aggregate array
        numObs = 0;
        for iTarget = 1:numTargets,
            for iZPos = 1:length(arrayZSlices{iPatient,iTarget}),
               numObs = numObs + length( arrayLayerPixThickness{iPatient, iTarget, iDist}{iZPos} );
            end
        end
        arrayTempCombAllObs = zeros(numObs,1,'double');
        
        %populate the combined/aggregate array
        numStart = 1;
        for iTarget = 1:numTargets,
            for iZPos = 1:length(arrayZSlices{iPatient,iTarget}),
               numObs = length( arrayLayerPixThickness{iPatient, iTarget, iDist}{iZPos} );
               numEnd = numStart+numObs-1;
               arrayTempCombAllObs(numStart:numEnd) = arrayLayerPixThickness{iPatient, iTarget, iDist}{iZPos};
               numStart = numEnd+1;
            end
        end
        
        %convert from pixel distance to absolute distance
        arrayCombAllObs{iPatient,iDist} = arrayTempCombAllObs*arrayPixRes(iTarget,numPatient);
        
    end
    
end
    
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Create the output figure
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
 
%initialise the figure window
figOut = figure;
set(figOut, 'Position', arrayOutputPlotScreenPos);
set(figOut, 'PaperUnits', 'centimeters', 'PaperSize', arrayOutputPlotPaperSize, 'PaperPosition', arrayOutputPlotPaperPosition );
hold on;

%plot each of the measured distances in its own figure sub-panel
for iDist = 1:numDistToMeasure,
    
    %initialise the sub-panel
    subplot(2,2,iDist);
    hold on;
    
        
    %move through each patient
    for iPatient = 1:numPatients,
        
        %perform kernel smoothing over the combined data for this patient
        % and obtain a probability density function
        [arrayProbDens, arrayDist] = ksdensity(arrayCombAllObs{iPatient,iDist});
                
        %plot the smoothed pdf in the color specified for the patient
        plot(arrayDist, arrayProbDens, '-', 'Color', arrayPatientColors{iPatient}, 'LineWidth', 2);
        
    end
    
    %set the plot origin at (0,0) while retaining maximum values/axis
    % scaling
    arrayXLim = get(gca, 'XLim');
    set(gca, 'XLim', [0 arrayXLim(2)]);
    arrayYLim = get(gca, 'YLim');
    set(gca, 'YLim', [0 arrayYLim(2)]);
    
    %label the axes
    xlabel('Absolute distance ({\mu}m)', 'FontSize', 14, 'FontName', 'Arial');
    ylabel('Probability density function', 'FontSize', 14, 'FontName', 'Arial');
    
    %specify the font size for labelling the axes tick marks
    set(gca, 'FontSize', 14, 'FontName', 'Arial');
    
    %label the sub-panel
    text(0.65*arrayXLim(2),0.85*arrayYLim(2),arrayDistStrings{iDist}, 'FontSize', 14, 'FontName', 'Arial', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    text(-0.3*arrayXLim(2), 1.05*arrayYLim(2), arraySubFigStrings{iDist}, 'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold');
   
end
hold off;
print(figOut, '-r300', '-dpng', [ stringOutputDir '\output_tissue_and_layer_thickness.png' ]);
close(figOut);
