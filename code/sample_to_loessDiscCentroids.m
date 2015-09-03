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
% Last Updated: 03/09/15
%
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Input Parameters
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
 
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
arrayTargetPaths = { ['RAF'];
                     ['RAF_ph'];
                     ['MEK'];
                     ['MEK_ph'];
                     ['ERK'];
                     ['ERK_ph'] };
numTargets = length(arrayTargetPaths);

%the sub-cellular localisations examined within this data set
arrayLocStrings = { [ 'C' ];
                    [ 'N' ] };
numLocalisations = length(arrayLocStrings);   

numSpatialBins = 28;
arrayNormDistRatio = [1 4 2];
numLoessWindowSize = 0.5;


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
%confirm the relative path occurs as expected
strProcImgDataPath = [ stringMountedFolder 'processed' strFoldSep ];
strRawImgDataPath = [ stringMountedFolder 'image' strFoldSep ];
if (exist(strProcImgDataPath, 'dir') == 7) && (exist(strProcImgDataPath, 'dir') == 7),
    
    arrayAllDataByLoc = cell(numLocalisations, 1);

    for iLoc = 1:numLocalisations,
        arrayAllDataByLoc{iLoc} = zeros(numSpatialBins, numTargets, numPatients, 'double');
        for iTarget = 1:numTargets,
            for iPatient = 1:numPatients,
                numPatient = arrayPatients(iPatient);

                %load the sample_loc.mat array
                stringSegImgPath = [ strProcImgDataPath arrayTargetPaths{iTarget} strFoldSep 'Pat_' num2str(iPatient) strFoldSep ];
                stringSmpLocPath = [ strProcImgDataPath arrayTargetPaths{iTarget} strFoldSep 'Pat_' num2str(iPatient) strFoldSep arrayLocStrings{iLoc} strFoldSep ];
                load([stringSmpLocPath 'sample_loc.mat']);
                structSampleLocs = SampleOutput;
                clear SampleOutput;

                %determine the path of the image stacks
                stringImgDataPath = [ strRawImgDataPath arrayTargetPaths{iTarget} strFoldSep 'Pat_' num2str(iPatient) strFoldSep 'image_data' strFoldSep ];

                %produce a 'SampleAnalysis' structured array
                structSignalIntensityData = produceSmpAna(structSampleLocs, arraySamplingKernel, stringImgDataPath, stringSegImgPath);

                numSamplePoints = length(structSignalIntensityData);
                arrayXToSmooth = zeros(numSamplePoints,1,'double');
                arrayYToSmooth = zeros(numSamplePoints,1,'uint8');
                for iSample = 1:numSamplePoints,
                    arrayXToSmooth(iSample) = structSignalIntensityData(iSample).NormDist;
                    arrayYToSmooth(iSample) = structSignalIntensityData(iSample).SigInt;
                end

                arrayLoessCurve = calculateLoessCurve( arrayXToSmooth, arrayYToSmooth, numSpatialBins, arrayNormDistRatio, numLoessWindowSize );

                arrayBinCentres = interp1(arrayLoessCurve(1,:), arrayLoessCurve(2,:), 0.5:1:(numSpatialBins-0.5));

                arrayAllDataByLoc{iLoc}(:,iTarget,iPatient) = arrayBinCentres;

            end
        end
    end
    
else
    disp('warning: data files cannot be located');
    
end
