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
stringDataPath = 'C:\wc\2015_epidermal_data\data';
strProcImgDataPath = [ stringDataPath strFoldSep 'processed' strFoldSep ];
strRawImgDataPath = [ stringDataPath strFoldSep 'image' strFoldSep ];

for iTarget = 1:numTargets,
    for iPatient = 1:numPatients,
        numPatient = arrayPatients(iPatient);

        for iLoc = 1:numLocalisations,
            %load the sample_loc.mat array
            stringSegImgPath = [ strProcImgDataPath arrayTargetPaths{iTarget} '\Pat_' num2str(iPatient) '\' ];
            stringSmpLocPath = [ strProcImgDataPath arrayTargetPaths{iTarget} '\Pat_' num2str(iPatient) '\' arrayLocStrings{iLoc} '\' ];
            load([stringSmpLocPath 'sample_loc.mat']);
            structSampleLocs = SampleOutput;
            clear SampleOutput;

            %determine the path of the image stacks
            stringImgDataPath = [ strRawImgDataPath arrayTargetPaths{iTarget} '\Pat_' num2str(iPatient) '\image_data\' ];

            %produce a 'SampleAnalysis' structured array
            structSignalIntensityData = produceSmpAna(structSampleLocs, arraySamplingKernel, stringImgDataPath, stringSegImgPath);

            [ arrayLoessCurve ] = calculateLoessCurve( arrayXToSmooth, arrayYToSmooth, numSpatialBins, arrayNormDistRatio, numLoessWindowSize )
            
            a=1;
        end
    end
end