%% sample_to_loessDiscCentroids.m
% This MATLAB script 
%
% For the normalized Hill differential equation model, parameterisation was
%  performed to optimise the model fit for the cytoplasmic and nuclear
%  phospho-ERK-1/2 data across all three patients; while a more 
%  comprehensive data set containing phospho-Raf-1 and phospho-MEK-1/2 was
%  compared against the resulting model fit. These are output as separate
%  files, created at the end of this script, which can then used by
%  subsequent MATLAB scripts.
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
% Last Updated: 11/11/15
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

%target proteins/phospho-proteins to be exported from the data
arrayTargetPaths = { ['CALM'];
                     ['ERK_ph'];
                     ['MEK_ph'];
                     ['RAF_ph'] };
numTargets = length(arrayTargetPaths);

%the sub-cellular localisations examined within this data set (note that
% when there are no sampled data present; e.g. for plasma-membrane
% phospho-MEK-1/2, an array of zeroes will be present)
arrayLocStrings = { [ 'C' ];
                    [ 'N' ];
                    [ 'M' ] };
numLocalisations = length(arrayLocStrings);   

numSpatialBins = 7;
arrayNormDistRatio = [ 1 4 2 ];
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
stringFormatForSampleLocationDataToLoad = 'MATLAB';

%sample analysis data
flagLoadSampleAnalysisData = true;
stringFormatForSampleAnalysisDataToLoad = 'MATLAB';
flagReCalcSampleAnalysisData = false;
flagSaveSampleAnalysisData = false;
stringFormatForSampleAnalysisDataToSave = 'MATLAB';

%loess curve data
flagLoadLoessCurveData = true;
stringFormatForLoessCurveDataToLoad = 'MATLAB';
flagReCalcLoessCurveData = false;
flagSaveLoessCurveData = false;
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

stringOutputFolder = 'c:\test\';

 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Perform Pre-Processing
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
 
arrayLoessDataByLoc = cell(numLocalisations, 1);
arrayDataMeanByLoc = cell(numLocalisations, 1);
arrayDataStdDevByLoc = cell(numLocalisations, 1);

    
if (exist(strProcImgDataPath, 'dir') == 7) && (exist(strProcImgDataPath, 'dir') == 7),

    for iLoc = 1:numLocalisations,
        %as noted above, the output data are taken at d_norm = [0:0.5:7]
        arrayLoessDataByLoc{iLoc} = zeros((numSpatialBins*2)+1, numTargets, numPatients, 'double');
        arrayDataMeanByLoc{iLoc} = zeros(numTargets, numPatients, 'double');
        arrayDataStdDevByLoc{iLoc} = zeros(numTargets, numPatients, 'double');
        
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
                            %save the new 'SampleAnalysis' structured array if
                            % specified
                            if flagSaveSampleAnalysisData,
                                flagResult = saveSampAnalysis(structSignalIntensityData, stringSmpLocPath, stringFormatForSampleAnalysisDataToSave);
                            end
                        end
                    end

                    %extract the required information from the structured
                    % array into arrays for loess smoothing
                    numSamplePoints = length(structSignalIntensityData);
                    arrayXToSmooth = zeros(numSamplePoints,1,'double');
                    arrayYToSmooth = zeros(numSamplePoints,1,'uint8');
                    for iSample = 1:numSamplePoints,
                        arrayXToSmooth(iSample) = structSignalIntensityData(iSample).NormDist;
                        arrayYToSmooth(iSample) = structSignalIntensityData(iSample).SigInt;
                    end
                    
                    arrayDataMeanByLoc{iLoc}(iTarget,iPatient) = mean(double(arrayYToSmooth));
                    arrayDataStdDevByLoc{iLoc}(iTarget,iPatient) = std(double(arrayYToSmooth));

                    %load the loess-smoothed data
                    if flagLoadLoessCurveData,
                        arrayLoessCurve = loadLoessCurve(stringSmpLocPath, stringFormatForLoessCurveDataToLoad);
                        if (arrayLoessCurve == -1),
                            disp(['warning: flagLoadLoessCurveData is set to true, but the loess data for ' arrayTargetPaths{iTarget} ' in the ' arrayLocStrings{iLoc} ...
                                   ' compartment cannot be found in the specified format (' stringFormatForLoessCurveDataToLoad '), ' ...
                                   'attempting to re-calculate these data and output will automatically be saved, however run time will be significantly increased' ]);
                            arrayLoessCurve = calculateLoessCurve( arrayXToSmooth, arrayYToSmooth, numSpatialBins, arrayNormDistRatio, numLoessWindowSize );
                            %automatically save these in the 'expected
                            % input' format
                            flagResult = saveLoessCurve(arrayLoessCurve, stringSmpLocPath, stringFormatForLoessCurveDataToLoad);
                        end
                    else
                        %re-perform loess smoothing if specified
                        if flagReCalcLoessCurveData,
                            arrayLoessCurve = calculateLoessCurve( arrayXToSmooth, arrayYToSmooth, numSpatialBins, arrayNormDistRatio, numLoessWindowSize );
                            if flagSaveLoessCurveData,
                                flagResult = saveLoessCurve(arrayLoessCurve, stringSmpLocPath, stringFormatForLoessCurveDataToSave);
                            end
                        end
                    end

                    %confirm the scaling of the normalised distance
                    % co-ordinate
                    numMaxNormDist = ceil(max(arrayLoessCurve(1,:)));
                    if ((numMaxNormDist == 7) || (numMaxNormDist == 14) || (numMaxNormDist == 21) || (numMaxNormDist == 28)),
                        numNormDistSpacer = double(numMaxNormDist)/14;
                    else
                        disp('warning: maximum normalised distance an unexpected value (not a factor of 7 = sum(arrayNormDistRatio))');
                        numMaxNormDist = 28;
                    end
                    arrayOutputNormDistPositions = [0:numNormDistSpacer:double(numMaxNormDist)];
                    
                    arrayOutputPoints = interp1(arrayLoessCurve(1,:), arrayLoessCurve(2,:), arrayOutputNormDistPositions);
                    arrayLoessDataByLoc{iLoc}(:,iTarget,iPatient) = arrayOutputPoints;
                    
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
%output the cytoplasmic and nuclear phospho-ERK data for fitting 
filePhosphoERKCytoForFitting = fopen([stringOutputFolder 'phosphoERK_cyto_toFit.csv' ], 'w+'); 

%print the header
fprintf(filePhosphoERKCytoForFitting, ['Patient']);
for iXPos = 1:((numSpatialBins*2)+1),
    fprintf(filePhosphoERKCytoForFitting, [',x=' num2str(iXPos)]);
end
fprintf(filePhosphoERKCytoForFitting, ['\n']);

%print the data by patient
for iPatient = 1:numPatients,
    %normalise the data around the mean/standard deviation (convert to
    % z-score); iLoc = 1 == C; iTarget = 2 == ERK_ph
    arrayOutputData = arrayLoessDataByLoc{1}(:,2,iPatient);
    arrayOutputDataNorm = (arrayOutputData - arrayDataMeanByLoc{1}(2,iPatient))/arrayDataStdDevByLoc{1}(2,iPatient);
    %output the normalised data
    fprintf(filePhosphoERKCytoForFitting, ['Pat' num2str(iPatient)]);
    for iXPos = 1:((numSpatialBins*2)+1),
        fprintf(filePhosphoERKCytoForFitting, [',' num2str(arrayOutputDataNorm(iXPos),'%f')]);
    end
    %line feed
    if iPatient < numPatients,
        fprintf(filePhosphoERKCytoForFitting, ['\n']);
    end
end

fclose(filePhosphoERKCytoForFitting);