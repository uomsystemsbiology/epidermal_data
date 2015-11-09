%% recreate_pMEK_heterogeneity_fig.m
% This MATLAB script reproduces Fig. AF4.1 from the GigaScience Data Note:
%   Cursons et al. (2015). Spatially-transformed fluorescence image data 
%    for ERK-MAPK and selected proteins within human epidermis.
%    GigaScience. Submitted Sept 2015.
%    doi: not-yet-known
%
% A number of functions are used by this script, some of which have
%  dependencies upon MATLAB Toolboxes:
%   - Image Processing Toolbox: this script directly calls the bwlabel, 
%                                   regionprops, and bwmorph functions
%                               the calculateDistancesToBoundary function also
%                                   calls bwlabel and regionprops
%   - Neural Network Toolbox: dist is called directly by this script
%   - Statistics and Machine Learning Toolbox: this script directly calls
%                                   the ksdensity function
%
% This script also uses the 'violin' function from the MATLAB File
%   Exchange: File ID: #45134
%       http://www.mathworks.com/matlabcentral/fileexchange/45134-violin-plot-based-on-kernel-density-estimation
%       violin.m was created by Holger Hoffman: hhoffmann@uni-bonn.de
%
% This script was created by Joe Cursons at the University of Melbourne
%   Systems Biology Laboratory:
%       joseph.cursons@unimelb.edu.au
%
% Note that if you are using a UNIX system, alpha transparency and lighting
%  have been disabled, as this can crash the Virtual Reference Environment
%  (due to OpenGL conflicts). If you are using UNIX without a VRE and would
%  like to produce shinier surfaces, you will need to modify the
%  conditional loops around lines 456-466 and 500-510.
%
% Last Updated: 03/11/15
%
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Input Parameters
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  

%specify the prefix string for the image data to be analysed
strImagePrefix = '6_Series004_z004_';
strTarget = 'MEK_ph';
numPatient = 2;

%specify the channels to be analysed in this analysis
% NB: image data have a _ch<#> prefix, this should be 0 for all of the
%  single-target image data
arrayChannelNums = 0;

%specify parameters for the spatial partitioning of the tissue
numSpatialBins = 7;
arrayNormDistRatio = [1 4 2];

%specify the number of 'spatial bins' to use for the histogram surface
% rendering
numSamplesPerSpatPartition = 50;
numSpatialBinsForHistSurf = sum(arrayNormDistRatio)*numSamplesPerSpatPartition;
%specify the number of 'signal bins' to use for the histogram surface
% rendering
numSpatHistSigBins = 64;
%specify surface interpolation settings (to improve visual display)
numSurfInterp = 4;
%specify 'spacing' parameters for adjacent histograms
numHistSurfSpacerMult = 3;
numHistSurfSpacerAdd = ceil(numHistSurfSpacerMult/2);

%specify output plotting settings
arrayOutputPlotPaperSize = [ 24, 30 ];
arrayOutputPlotPaperPosition = [ 0, 0, 18, 24 ];
arrayOutputPlotScreenPos = [ 50 50 2450 3050 ];
numPlotFontSize = 12;

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
addpath(genpath(strCurrDir));

%manipulate the file path to determine the appropriate folders
arrayCurrDirFoldSepPos = strfind(strCurrDir, strFoldSep);

if arrayCurrDirFoldSepPos(end) == length(strCurrDir),
    %there is a backslash at the end
    strBaseDir = strCurrDir(1:(arrayCurrDirFoldSepPos(end-1)));
else
    strBaseDir = strCurrDir(1:(arrayCurrDirFoldSepPos(end)));
end

%set the relative image data path
stringImageDataPath = [ strBaseDir 'image' strFoldSep strTarget strFoldSep 'Pat_' num2str(numPatient) strFoldSep 'image_data' strFoldSep ];

%set the relative processed data path
stringProcessedDataPath = [ strBaseDir 'processed' strFoldSep strTarget strFoldSep 'Pat_' num2str(numPatient) strFoldSep  ];

%and just output to the root directory
stringOutputDir = strBaseDir;

 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Match the segmented objects and extract the pixel intensity data
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  

%determine the total number of channels
numChannels = length(arrayChannelNums);

%load the image data into an array for subsequent analysis
arrayInputImages = cell(numChannels,1);
for iChan = 1:numChannels,
    numChan = arrayChannelNums(iChan);
    strImageName = [ strImagePrefix 'ch' num2str(numChan, '%02i') '.tif' ];
    arrayInputImages{iChan} = imread([stringImageDataPath strImageName]);
end

%load the cytoplasmic and nuclear binarised mask files
imgBinMskCyto = imread([ stringProcessedDataPath strImagePrefix 'cyto.tif' ]);
imgBinMskCyto = imgBinMskCyto(:,:,1); %GIMP TIFF files have a transparency channel
imgBinMskNuc = imread([ stringProcessedDataPath strImagePrefix 'nuc.tif' ]);
imgBinMskNuc = imgBinMskNuc(:,:,1);

%to match the cytoplasmic and nuclear objects belonging to the same cell,
% I first use bwlabel to identify demarcated object (by '8' separation;
% i.e. no touching pixels)
imgLabCyto = bwlabel(imgBinMskCyto, 8);
imgLabNuc = bwlabel(imgBinMskNuc, 8);

%check the number of objects is equal across both masks
if ~(max(imgLabCyto(:)) == max(imgLabNuc(:))),
    disp('warning: the number of cytoplasmic and nuclear objects does not match');
end

%use the regionprops function to extract the centroids and pixel indices
% (as lists) of the labelled objects
structCytoProps = regionprops(imgLabCyto, 'Centroid', 'PixelIdxList');
structNucProps = regionprops(imgLabNuc, 'Centroid', 'PixelIdxList');

%confirm that the resulting data structures are also the same length
if ~(length(structCytoProps) == length(structNucProps)),
    disp('warning: the number of cytoplasmic and nuclear objects does not match');
end

%extract the centroids of cytoplasmic objects into an array
arrayCytoCentroids = zeros(2, length(structCytoProps), 'int16');
for iObj = 1:length(structCytoProps),
    arrayTempCent = structCytoProps(iObj).Centroid;
    arrayCytoCentroids(:,iObj) = arrayTempCent;
end

%extract the centroids of nuclear objects into an array
arrayNucCentroids = zeros(2, length(structNucProps), 'int16');
for iObj = 1:length(structNucProps),
    arrayTempCent = structNucProps(iObj).Centroid;
    arrayNucCentroids(:,iObj) = arrayTempCent;
end


%extract the signal intensity data from the cytoplasmic objects into arrays
arrayCytoSigInt = cell(length(structCytoProps),3);
arrayNucSigInt = cell(length(structNucProps),3);
%and create image files where cytoplasmic and nuclear objects are labelled
% by their 'object number' for curation of the cytoplasm:nucleus matching
imageLabelledCytoplasm = zeros(size(imgLabCyto,1), size(imgLabCyto,2), 'int16');
imageLabelledNuclei = zeros(size(imgLabNuc,1), size(imgLabNuc,2), 'int16');
%move through every nuclear object
for iObj = 1:length(structNucProps),
    
    %extract pixel co-ordinates of the object
    arrayTempPixList = structNucProps(iObj).PixelIdxList;
    imageLabelledNuclei(arrayTempPixList) = int16(iObj);
    for iChan = 1:numChannels,
        arrayNucSigInt{iObj,iChan} = arrayInputImages{iChan}(arrayTempPixList);
    end
    
    %extract centroid co-ordinates of the nucleus object into an array
    arrayTempCentroid = structNucProps(iObj).Centroid;
    
    %append the centroid co-ordinates of all the cytoplasmic objects
    arrayForDist = [arrayTempCentroid' arrayCytoCentroids];
    
    %calculate distances between  pixel co-ordinates of the object using
    % the dist function 
    arrayDist = dist(arrayForDist);
    %identify the closest cytoplasmic object
    [~,numClosestCyto] = min(arrayDist(1, 2:end));
    
    %extract data from the 'matched' cytoplasmic object
    arrayTempPixList = structCytoProps(numClosestCyto).PixelIdxList;
    imageLabelledCytoplasm(arrayTempPixList) = int16(iObj);
    for iChan = 1:numChannels,
        arrayCytoSigInt{iObj,iChan} = arrayInputImages{iChan}(arrayTempPixList);
    end
end

%determine the final number of output objects
numOutputObjects = max(imageLabelledNuclei(:));

%create an output array for the normalised-distance value for each of the
% output nuclei positions
arrayNormDist = zeros(numOutputObjects,1,'double');

%load the corresponding segmentation data for the tissue
imageTissueLayers = imread([ stringProcessedDataPath 'im_layers_z5.tiff' ]);

%move through each of the output objects and calculate the normalised
%distance co-ordinate
for iObj = 1:numOutputObjects,
    %figure out where the sample lies within the tissue
    if (imageTissueLayers(int16(structNucProps(iObj).Centroid(2)), int16(structNucProps(iObj).Centroid(1))) == 1),

        %extract the corresponding tissue layer boundaries
        imageBasalLamina = imread([ stringProcessedDataPath 'im_basal_memb_z5.tiff' ]);
        imageBoundOne = imread([ stringProcessedDataPath 'im_bound1_z5.tiff' ]);
        %calculate the cellular distances
        numBasalLaminaDistance = calculateDistancesToBoundary(imageBasalLamina, structNucProps(iObj).Centroid');
        numBoundaryOneDistance = calculateDistancesToBoundary(imageBoundOne, structNucProps(iObj).Centroid');
        %add the appropriate offset and perform linear interpolation
        arrayNormDist(iObj) = 0 + (  numBasalLaminaDistance/(numBasalLaminaDistance + numBoundaryOneDistance)  );

    elseif (imageTissueLayers(int16(structNucProps(iObj).Centroid(2)), int16(structNucProps(iObj).Centroid(1))) == 2),

        %extract the corresponding tissue layer boundaries
        imageBoundOne = imread([ stringProcessedDataPath 'im_bound1_z5.tiff' ]);
        imageBoundTwo = imread([ stringProcessedDataPath 'im_bound2_z5.tiff' ]);
        %calculate the cellular distances
        numBoundaryOneDistance = calculateDistancesToBoundary(imageBoundOne, structNucProps(iObj).Centroid');
        numBoundaryTwoDistance = calculateDistancesToBoundary(imageBoundTwo, structNucProps(iObj).Centroid');
        %add the appropriate offset and perform linear interpolation
        arrayNormDist(iObj) = 1 + (  numBoundaryOneDistance/(numBoundaryOneDistance + numBoundaryTwoDistance)  );

    elseif (imageTissueLayers(int16(structNucProps(iObj).Centroid(2)), int16(structNucProps(iObj).Centroid(1))) == 3),

        %extract the corresponding tissue layer boundaries
        imageBoundTwo = imread([ stringProcessedDataPath 'im_bound2_z5.tiff' ]);
        imageOuterBound = imread([ stringProcessedDataPath 'im_outer_bound_z5.tiff' ]);
        %calculate the cellular distances
        numBoundaryTwoDistance = calculateDistancesToBoundary(imageBoundTwo, structNucProps(iObj).Centroid');
        numOuterBoundaryDistance = calculateDistancesToBoundary(imageOuterBound, structNucProps(iObj).Centroid');
        %add the appropriate offset and perform linear interpolation
        arrayNormDist(iObj) = 2 + (  numBoundaryTwoDistance/(numBoundaryTwoDistance + numOuterBoundaryDistance)  );

    %and samples outside the tissue
    elseif (imageTissueLayers(int16(structNucProps(iObj).Centroid(2)), int16(structNucProps(iObj).Centroid(1))) == -1),
        disp('warning: nucleus centroid appears to be located within the dermis');
    elseif (imageTissueLayers(int16(structNucProps(iObj).Centroid(2)), int16(structNucProps(iObj).Centroid(2))) == -3),
        disp('warning: nucleus centroid appears to be located outside of the tissue');
    else
        disp('warning: the layer image may not be loading properly');
    end
end


 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Produce output plot to curate the object alignment (matching nucleus 
%   and cytoplasm masks)
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%create the color map depending on the number of objects to be plotted
arrayColorMap = jet(double(numOutputObjects));

%create the output image as an RGB, 8-bit array
imageTestSeg = zeros(size(imgLabCyto,1), size(imgLabCyto,2),3, 'uint8');
for iObj = 1:numOutputObjects,
    %identify the pixel co-ordinates of the object
    arrayTempNucObj = (imageLabelledNuclei == iObj);
    arrayTempCytoObj = (imageLabelledCytoplasm == iObj);
    
    %use the bwmorph function to remove 'internal pixels'
    imageNucObjOutline = bwmorph(arrayTempNucObj, 'remove');
    imageCytoObjOutline = bwmorph(arrayTempCytoObj, 'remove');
    
    %combine the nuclear and cytoplasmic objects
    imageCombOutline = imageNucObjOutline | imageCytoObjOutline;
    
    %extract pixel-coordinates of the combined cyto/nuc outlines
    structOutlineProps = regionprops(imageCombOutline, 'PixelList');
    %draw their outlines on the image, using the specified colour map
    for jObj = 1:length(structOutlineProps),
        for iPix = 1:length(structOutlineProps(jObj).PixelList),
            numX = structOutlineProps(jObj).PixelList(iPix,1);
            numY = structOutlineProps(jObj).PixelList(iPix,2);
            imageTestSeg(numY, numX,:) = uint8(arrayColorMap(iObj,:)*255);
        end
    end
end

%create the output image
imwrite(imageTestSeg, [stringOutputDir 'outputCurateCytoNucAlignment.tiff' ]);

 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Produces 'histogram surface' for the output plots
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  

%calculate spatial discretisation parameters from the input settings
% defined above
[ arrayIntervalCentres, arrayDistIntervals, arrayDivisionIndices ] = calculateSpatialDivisions( numSpatialBins, arrayNormDistRatio );
arrayLayerBoundaries = arrayDivisionIndices-1;

%create parameters required for output which are dependent upon the data
numMaxCytoNodeSignal = 0;
numMaxNucNodeSignal = 0;

%construct the histogram arrays for the data
arrayCytoHistSurf = zeros(numSpatHistSigBins, numSpatialBinsForHistSurf*numHistSurfSpacerMult);
arrayNucHistSurf = zeros(numSpatHistSigBins, numSpatialBinsForHistSurf*numHistSurfSpacerMult);
for iObj = 1:numOutputObjects,

    %adjust the max signal parameters if necessary
    numMaxCytoNodeSignal = max([arrayCytoSigInt{iObj,1};numMaxCytoNodeSignal]);
    numMaxNucNodeSignal = max([arrayNucSigInt{iObj,1};numMaxNucNodeSignal]);
    
    %check the normalised-distance position of the object
    if (arrayNormDist(iObj) >= 0) && (arrayNormDist(iObj) < 1),
        numSpatBin = int64(arrayNormDistRatio(1)*numSamplesPerSpatPartition*arrayNormDist(iObj));
    elseif (arrayNormDist(iObj) >= 1) && (arrayNormDist(iObj) < 2),
        numSpatBin = int64(arrayNormDistRatio(2)*numSamplesPerSpatPartition*(arrayNormDist(iObj)-1) + arrayNormDistRatio(1)*numSamplesPerSpatPartition);
    elseif (arrayNormDist(iObj) >= 2) && (arrayNormDist(iObj) <= 3),
        numSpatBin = int64(arrayNormDistRatio(3)*numSamplesPerSpatPartition*(arrayNormDist(iObj)-2) + sum(arrayNormDistRatio(1:2))*numSamplesPerSpatPartition);
    else
        disp('warning: object has a normalised distance value less than zero or greater than three');
    end

    %examine the pixel intensities to construct the histograms
    arrayBinEdges = -0.5:256/numSpatHistSigBins:255.5;
    arrayBinCent = arrayBinEdges(1:(end-1))+1.5;
    arrayCytoBinFreq = histc(arrayCytoSigInt{iObj,1}, arrayBinEdges);
    arrayCytoBinFreq = arrayCytoBinFreq(1:(end-1)); %remove the final bin (freq > highest edge)
    arrayNucBinFreq = histc(arrayNucSigInt{iObj,1}, arrayBinEdges);
    arrayNucBinFreq = arrayNucBinFreq(1:(end-1)); %remove the final bin (freq > highest edge)

    %convert into a spaced surface to render
    arrayCytoHistSurf(1:numSpatHistSigBins, numSpatBin*numHistSurfSpacerMult+numHistSurfSpacerAdd) = arrayCytoBinFreq;
    arrayNucHistSurf(1:numSpatHistSigBins, numSpatBin*numHistSurfSpacerMult+numHistSurfSpacerAdd) = arrayNucBinFreq;
end
%interpolate the surfaces
arrayCytoHistSurf = interp2(arrayCytoHistSurf);
arrayNucHistSurf = interp2(arrayNucHistSurf);
%calculate the alpha (transparency) as a log transform (so very low values
% are 'see through')
arrayCytoHistSurfAlpha = log(arrayCytoHistSurf);
arrayNucHistSurfAlpha = log(arrayNucHistSurf);

%sort the objects by their normalised distance co-ordinates
[~,arraySortByNormDistIndex] = sort(arrayNormDist);

%extract these data into arrays to create the violin plots
arrayCytoForViolinPlot = cell(1,numOutputObjects);
arrayNucForViolinPlot = cell(1,numOutputObjects);
arrayXLabs = cell(1,numOutputObjects);
for iObj = 1:numOutputObjects,
    arrayCytoForViolinPlot{iObj} = zeros(length(arrayCytoSigInt{arraySortByNormDistIndex(iObj),1}),1);
    arrayCytoForViolinPlot{iObj}(:) = arrayCytoSigInt{arraySortByNormDistIndex(iObj),1};
    arrayNucForViolinPlot{iObj} = zeros(length(arrayNucSigInt{arraySortByNormDistIndex(iObj),1}),1);
    arrayNucForViolinPlot{iObj}(:) = arrayNucSigInt{arraySortByNormDistIndex(iObj),1};
    arrayXLabs{iObj} = ['d_norm=' num2str(arrayNormDist(arraySortByNormDistIndex(iObj)))];
end


 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Produce the Output Figure
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%create the figure and specify its size/position
figOut = figure;
set(figOut, 'Position', arrayOutputPlotScreenPos);
set(figOut, 'PaperUnits', 'centimeters', 'PaperSize', arrayOutputPlotPaperSize, 'PaperPosition', arrayOutputPlotPaperPosition );

%create the cytoplasmic ksdensity plots which span panels 1 & 3 of a 
% 5*2 subplot array
subplot(10,2,[1,3]);
hold on;
%for each output object
for iObj = 1:numOutputObjects,
    
    %as sorted by normalised distance
    numObj = arraySortByNormDistIndex(iObj);
    
    %determine the probability density
    [arrayProbDens, arrayDist] = ksdensity(double(arrayCytoSigInt{numObj,1}));
    
    %and output in the corresponding colour
    plot(arrayDist, arrayProbDens, '-', 'Color', uint8(arrayColorMap(numObj,:)*255));
            
end
hold off;
%label axes and specify ranges
set(gca, 'XLim', [0 255]);
xlabel('Cytoplasmic Signal Intensity');
ylabel('Prob. Dens.');
%shift the plot slightly to improve overall figure layout..
arrayCytoPDFPlotPos = get(gca, 'Position');
set(gca, 'Position', (arrayCytoPDFPlotPos + [ 0.00 0.05 0.00 0.00 ]));

%create the cytoplasmic ksdensity plots which span panels 2 & 4 of a 
% 5*2 subplot array
subplot(10,2,[2,4]);
hold on;
%for each output object
for iObj = 1:numOutputObjects,
    
    %as sorted by normalised distance
    numObj = arraySortByNormDistIndex(iObj);
    
    %determine the probability density
    [arrayProbDens, arrayDist] = ksdensity(double(arrayNucSigInt{numObj,1}));

    %and output in the corresponding colour
    plot(arrayDist, arrayProbDens, '-', 'Color', uint8(arrayColorMap(numObj,:)*255));
            
end
hold off;
%label axes and specify ranges
set(gca, 'XLim', [0 255]);
xlabel('Nuclear Signal Intensity');
ylabel('Prob. Dens.');
%shift the plot slightly to improve overall figure layout..
arrayNucPDFPlotPos = get(gca, 'Position');
set(gca, 'Position', (arrayNucPDFPlotPos + [ 0.00 0.05 0.00 0.00 ]));


%create the cytoplasmic violin plots which span panels 5-8 of a 
% 5*2 subplot array
subplot(10,2,5:8);
%create the violin plots of the cytoplasmic data
[arrayHandles,L,MX,MED,bw]=violin(arrayCytoForViolinPlot);
%label the axes
xlabel('cells ordered by d_{norm}');
ylabel('Cyto. Sig. Int.');

%create the nuclear violin plots which span panels 9-12 of a 5*2 subplot 
% array
subplot(10,2,9:12);
%create the violin plots of the cytoplasmic data
[arrayHandles,L,MX,MED,bw]=violin(arrayNucForViolinPlot);
%label the axes
xlabel('cells ordered by d_{norm}');
ylabel('Nuc. Sig. Int.');


%create the cytoplasmic intensity histogram surface plots which span panels 
% 13-16 of a 5*2 subplot array
subplot(10,2,13:16);
hold on;
%plot the cytoplasmic data spatial histograms
if isdeployed,
    %don't use alpha data or lighting, because this can crash the virtual
    % box (due to OpenGL rendering issues)
    surf(arrayCytoHistSurf, 'FaceColor','interp',  'EdgeColor','none');
else
    %plot the surface using alpha transparency
    surf(arrayCytoHistSurf, ...
        'FaceColor','interp',  'EdgeColor','none',  'FaceAlpha','interp',...
        'AlphaDataMapping','scaled',  'AlphaData',arrayCytoHistSurfAlpha, ...
        'FaceLighting','phong', 'AmbientStrength',0.8);
end
%format the plot
axis([0, size(arrayCytoHistSurf,2), 0, size(arrayCytoHistSurf,1), 0, max(arrayCytoHistSurf(:))*1.1]);
view([15 64]);
set(gca, 'FontName', 'Times New Roman', 'FontSize',  numPlotFontSize);
text(0.5*double(numSpatialBinsForHistSurf*2*numHistSurfSpacerMult), -0.3*numSpatHistSigBins, 'Normalized Distance (d_{norm})', 'FontName', 'Times New Roman', 'FontSize',  numPlotFontSize, 'HorizontalAlignment', 'center');
for iSpatPart = 1:sum(arrayNormDistRatio),
    text(double(iSpatPart-0.5)*numSpatialBinsForHistSurf*2*numHistSurfSpacerMult, -2, 0, num2str(iSpatPart))
end
set(gca, 'XTick', ((arrayLayerBoundaries*numSamplesPerSpatPartition*2*numHistSurfSpacerMult)+1));
set(gca, 'XTickLabel', []);
text(1.2*double(numSpatialBinsForHistSurf*2*numHistSurfSpacerMult), numSpatHistSigBins*0.4, {'signal';'intensity'}, 'FontName', 'Times New Roman', 'FontSize',  numPlotFontSize, 'HorizontalAlignment', 'center');
set(gca, 'YTick', [ 0 numMaxCytoNodeSignal ]);
set(gca, 'YTickLabel', []);
text(-0.5*double(arrayLayerBoundaries(2)*2*numHistSurfSpacerMult), 0, 0, '0', 'FontName', 'Times New Roman', 'FontSize',  numPlotFontSize);
text(1.03*double(numSpatialBinsForHistSurf*numHistSurfSpacerMult*2), numSpatHistSigBins*2, 0, num2str(numMaxCytoNodeSignal-1), 'FontName', 'Times New Roman', 'FontSize',  numPlotFontSize);
text(1.03*double(numSpatialBinsForHistSurf*numHistSurfSpacerMult*2), 0.5, 0, '0', 'FontName', 'Times New Roman', 'FontSize',  numPlotFontSize);
text(-0.25*double(numSpatialBinsForHistSurf*numHistSurfSpacerMult),0, 0.5*double(max(arrayCytoHistSurf(:))), 'freq.', 'FontName', 'Times New Roman', 'FontSize',  numPlotFontSize);
set(gca, 'FontName', 'Times New Roman', 'FontSize',  numPlotFontSize);
numRoundedZTick = double(int16(max(arrayCytoHistSurf(:))/100)*100);
set(gca, 'ZTick', [0 numRoundedZTick] );
set(gca, 'ZTickLabel', [] );
text(-0.17*double(numSpatialBinsForHistSurf*numHistSurfSpacerMult), -0.01*numSpatHistSigBins, numRoundedZTick, num2str(numRoundedZTick), 'FontName', 'Times New Roman', 'FontSize',  numPlotFontSize, 'VerticalAlignment', 'middle');

hold off;
handleCytoPDFSurfPlot = gca;
arrayCytoPDFSurfPlotPos = get(gca, 'Position');


%create the nuclear intensity histogram surface plots which span panels 
% 17-20 of a 5*2 subplot array
subplot(10,2,17:20);
hold on;
%plot the nuclear data spatial histograms
if isdeployed,
    %don't use alpha data or lighting, because this can crash the virtual
    % box (due to OpenGL rendering issues)
    surf(arrayNucHistSurf, 'FaceColor','interp',  'EdgeColor','none');
else
    %plot the surface using alpha transparency
    surf(arrayNucHistSurf, ...
        'FaceColor','interp',  'EdgeColor','none',  'FaceAlpha','interp',...
        'AlphaDataMapping','scaled',  'AlphaData',arrayNucHistSurfAlpha, ...
        'FaceLighting','phong', 'AmbientStrength',0.8);
end
%format the plot
axis([0, size(arrayNucHistSurf,2), 0, size(arrayNucHistSurf,1), 0, max(arrayNucHistSurf(:))*1.1]);
view([15 64]);
set(gca, 'FontName', 'Times New Roman', 'FontSize',  numPlotFontSize);
text(0.5*double(numSpatialBinsForHistSurf*2*numHistSurfSpacerMult), -0.3*numSpatHistSigBins, 'Normalized Distance (d_{norm})', 'FontName', 'Times New Roman', 'FontSize',  numPlotFontSize, 'HorizontalAlignment', 'center');
for iSpatPart = 1:sum(arrayNormDistRatio),
    text(double(iSpatPart-0.5)*numSpatialBinsForHistSurf*2*numHistSurfSpacerMult, -2, 0, num2str(iSpatPart))
end
set(gca, 'XTick', ((arrayLayerBoundaries*numSamplesPerSpatPartition*2*numHistSurfSpacerMult)+1));
set(gca, 'XTickLabel', []);
text(1.2*double(numSpatialBinsForHistSurf*2*numHistSurfSpacerMult), numSpatHistSigBins*0.4, {'signal';'intensity'}, 'FontName', 'Times New Roman', 'FontSize',  numPlotFontSize, 'HorizontalAlignment', 'center');
set(gca, 'YTick', [ 0 numMaxNucNodeSignal ]);
set(gca, 'YTickLabel', []);
text(-0.5*double(arrayLayerBoundaries(2)*2*numHistSurfSpacerMult), 0, 0, '0', 'FontName', 'Times New Roman', 'FontSize',  numPlotFontSize);
text(1.03*double(numSpatialBinsForHistSurf*numHistSurfSpacerMult*2), numSpatHistSigBins*2, 0, num2str(numMaxNucNodeSignal-1), 'FontName', 'Times New Roman', 'FontSize',  numPlotFontSize);
text(1.03*double(numSpatialBinsForHistSurf*numHistSurfSpacerMult*2), 0.5, 0, '0', 'FontName', 'Times New Roman', 'FontSize',  numPlotFontSize);
text(-0.25*double(numSpatialBinsForHistSurf*numHistSurfSpacerMult),0, 0.5*double(max(arrayNucHistSurf(:))), 'freq.', 'FontName', 'Times New Roman', 'FontSize',  numPlotFontSize);
set(gca, 'FontName', 'Times New Roman', 'FontSize',  numPlotFontSize);
numRoundedZTick = double(int16(max(arrayNucHistSurf(:))/100)*100);
set(gca, 'ZTick', [0 numRoundedZTick] );
set(gca, 'ZTickLabel', [] );
text(-0.17*double(numSpatialBinsForHistSurf*numHistSurfSpacerMult), -0.01*numSpatHistSigBins, numRoundedZTick, num2str(numRoundedZTick), 'FontName', 'Times New Roman', 'FontSize',  numPlotFontSize, 'VerticalAlignment', 'middle');

%move the nuclear plot slightly to improve appearance
arrayNucPDFSurfPlotPos = get(gca, 'Position');
set(gca, 'Position', (arrayNucPDFSurfPlotPos + [ 0.00 -0.05 0.00 0.00 ]));

hold off;

%move the cytoplasmic plot slightly to improve appearance
% NB: because this moves the plot down, it needs to be moved after the
% nuclear plot, or it will disappear when the nuclear histogram surface is
% rendered
set(handleCytoPDFSurfPlot, 'Position', (arrayCytoPDFSurfPlotPos + [ 0.00 -0.05 0.00 0.00 ]));

%print to the specified output file
print(figOut, '-r300', '-dpng', [ stringOutputDir 'FullCellSeg_pMEK_IntDistPDFs.png' ]);

%and close the figure window
close(figOut);