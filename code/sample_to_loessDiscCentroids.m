%% sample_to_loessDiscCentroids.m
% This MATLAB script reads in the sample location data, extracts pixel
%  intensities from the raw image data, maps to the normalised distance
%  co-ordinate, performs loess smoothing, and outputs .csv/.tiff at various
%  stages.
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
% This script was created by Joe Cursons at the University of Melbourne
%   Systems Biology Laboratory:
%       joseph.cursons@unimelb.edu.au
%
% Last Updated: 22/09/15
% 
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Input Parameters
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%define the local neighbourhood of pixels to sample the signal intensity 
% for each region-of-interest
arraySamplingKernel = [ 1 1 1 1 1;      
                        1 1 1 1 1;
                        1 1 1 1 1;
                        1 1 1 1 1;
                        1 1 1 1 1 ];

%the number of patients contained within this data set
arrayPatients = 1:3;
numPatients = length(arrayPatients);


%the target proteins/phospho-proteins contained within this data set
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
% arrayTargetPaths = { ['CALM'];
%                      ['ERK_ph'];
%                      ['MEK_ph'];
%                      ['RAF_ph'] };
numTargets = length(arrayTargetPaths);

%the sub-cellular localisations examined within this data set
arrayLocStrings = { [ 'C' ];
                    [ 'N' ];
                    [ 'M' ] };
numLocalisations = length(arrayLocStrings);   

numSpatialBins = 28;
arrayNormDistRatio = [1 4 2];
numLoessWindowSize = 0.5;


 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Specify Flags/Strings to Control Data Input/Output and Formats
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%different stages of the data processing script have different data formats
% available (specified below), but a list of strings which will be accepted
% by the data processing functions include:
%   MATLAB data files (.mat):
%       'MATLAB', '.mat', 'matlab'
%   Image files (.tiff/.TIFF):
%       'TIFF', '.TIFF', 'tif', '.tif', 'tiff', '.tiff', 'image', 'Image'
 
%sample location data - MATLAB or 
flagLoadSampleLocationData = true;
stringFormatForSampleLocationDataToLoad = 'Image';

%sample analysis data
flagLoadSampleAnalysisData = true;
stringFormatForSampleAnalysisDataToLoad = 'MATLAB';
flagReCalcSampleAnalysisData = false;
flagSaveSampleAnalysisData = false;
stringFormatForSampleAnalysisDataToSave = 'MATLAB';

%loess curve data
flagLoadLoessCurveData = false;
stringFormatForLoessCurveDataToLoad = 'MATLAB';
flagReCalcLoessCurveData = true;
flagSaveLoessCurveData = true;
stringFormatForLoessCurveDataToSave = 'MATLAB';

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

%load shared_functions
addpath(genpath([strCurrDir strFoldSep 'functions']));

%and the relative path of the data
arrayCurrDirSeps = strfind(strCurrDir, strFoldSep);
numDirSeps = length(arrayCurrDirSeps);
stringMountedFolder = strCurrDir(1:arrayCurrDirSeps(numDirSeps-1));
stringMountedFolder = 'C:\wc\2015_epidermal_data\data\';
%confirm the relative path occurs as expected
strProcImgDataPath = [ stringMountedFolder 'processed' strFoldSep ];
strRawImgDataPath = [ stringMountedFolder 'image' strFoldSep ];


 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Perform Pre-Processing
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
 
if (exist(strProcImgDataPath, 'dir') == 7) && (exist(strProcImgDataPath, 'dir') == 7),
    
    arrayAllDataByLoc = cell(numLocalisations, 1);

    for iLoc = 1:numLocalisations,
        arrayAllDataByLoc{iLoc} = zeros(numSpatialBins, numTargets, numPatients, 'double');
        for iTarget = 1:numTargets,
            for iPatient = 1:numPatients,
                numPatient = arrayPatients(iPatient);

                stringSegImgPath = [ strProcImgDataPath arrayTargetPaths{iTarget} strFoldSep 'Pat_' num2str(iPatient) strFoldSep ];
                stringSmpLocPath = [ strProcImgDataPath arrayTargetPaths{iTarget} strFoldSep 'Pat_' num2str(iPatient) strFoldSep arrayLocStrings{iLoc} strFoldSep ];
                
                %check that the specified sub-cellular localisation exists
                % for data processing
                if (exist(stringSmpLocPath, 'dir') == 7) && (exist(stringSmpLocPath, 'dir') == 7),

                    %load the sample location data
                    if flagLoadSampleLocationData,
                        structSampleLocs = loadSampLocs(stringSmpLocPath, stringFormatForSampleLocationDataToLoad);
%                         structSampleLocsFromMATLAB = loadSampLocs(stringSmpLocPath, 'MATLAB');
                    else
                        %there are some issues with including the scripts for
                        % performing new sampling as it requires the MATLAB GUI
                        % and this cannot be distributed with the VRE (which
                        % uses the MATLAB component run time). JC is working
                        % on an alternative.
                        disp('warning: code has not yet been integrated for writing new sample location data');
                    end

                    %determine the path of the image stacks
                    stringImgDataPath = [ strRawImgDataPath arrayTargetPaths{iTarget} strFoldSep 'Pat_' num2str(iPatient) strFoldSep 'image_data' strFoldSep ];

                    %load the sample analysis data
                    if flagLoadSampleAnalysisData,
                        structSignalIntensityData = loadSampAnalysis(stringSmpLocPath, stringFormatForSampleAnalysisDataToLoad);
                    else
                        %produce a 'SampleAnalysis' structured array
                        if flagReCalcSampleAnalysisData,
                            structSignalIntensityData = produceSmpAna(structSampleLocs, arraySamplingKernel, stringImgDataPath, stringSegImgPath);
%                             structSignalIntensityDataFromMATLAB = produceSmpAna(structSampleLocsFromMATLAB, arraySamplingKernel, stringImgDataPath, stringSegImgPath);
                            %save the new 'SampleAnalysis' structured array if
                            % specified
                            if flagSaveSampleAnalysisData,
                                flagResult = saveSampAnalysis(structSignalIntensityData, stringSmpLocPath, stringFormatForSampleAnalysisDataToSave);
                            end
                        end
                    end

                    numSamplePoints = length(structSignalIntensityData);
                    arrayXToSmooth = zeros(numSamplePoints,1,'double');
%                     arrayYToSmooth = zeros(numSamplePoints,1,'uint8');
%                     arrayXToSmoothFromMATLAB = zeros(numSamplePoints,1,'double');
%                     arrayYToSmoothFromMATLAB = zeros(numSamplePoints,1,'uint8');
                    for iSample = 1:numSamplePoints,
                        arrayXToSmooth(iSample) = structSignalIntensityData(iSample).NormDist;
                        arrayYToSmooth(iSample) = structSignalIntensityData(iSample).SigInt;
%                         arrayXToSmoothFromMATLAB(iSample) = structSignalIntensityDataFromMATLAB(iSample).NormDist;
%                         arrayYToSmoothFromMATLAB(iSample) = structSignalIntensityDataFromMATLAB(iSample).SigInt;
                    end

                    %load the loess-smoothed data
                    if flagLoadLoessCurveData,
                        arrayLoessCurve = loadLoessCurve(stringSmpLocPath, stringFormatForLoessCurveDataToLoad);
                    else
                        %re-perform loess smoothing if specified
                        if flagReCalcLoessCurveData,
                            arrayLoessCurve = calculateLoessCurve( arrayXToSmooth, arrayYToSmooth, numSpatialBins, arrayNormDistRatio, numLoessWindowSize );
%                             arrayLoessCurveFromMATLAB = calculateLoessCurve( arrayXToSmoothFromMATLAB, arrayYToSmoothFromMATLAB, numSpatialBins, arrayNormDistRatio, numLoessWindowSize );
                            if flagSaveLoessCurveData,
                                flagResult = saveLoessCurve(arrayLoessCurve, stringSmpLocPath, stringFormatForLoessCurveDataToSave);
                            end
                        end
                    end

                    arrayBinCentres = interp1(arrayLoessCurve(1,:), arrayLoessCurve(2,:), 0.5:1:(numSpatialBins-0.5));

                    arrayAllDataByLoc{iLoc}(:,iTarget,iPatient) = arrayBinCentres;
                else
                    disp(['warning: ' arrayTargetPaths{iTarget} ', Patient ' num2str(iPatient) ' has no data for localisation ' arrayLocStrings{iLoc}]);
                end

            end
        end
    end
    
else
    disp('warning: data files cannot be located');
    
end


 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Create Output
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  