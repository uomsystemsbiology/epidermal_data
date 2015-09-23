%% sample_to_loessDiscCentroids.m
% This MATLAB script moves through the complete set of processed data, 
%  extracts sample location data from the MATLAB data format (.mat), and 
%  then saves them using a TIFF format which is more amenable to analyses 
%  with open source languages such as python.
%
% A number of functions are used by this script, some of which have
%  dependencies upon MATLAB Toolboxes:
%   - Image Processing Toolbox: this script calls loadSampLocs and
%                                   saveSampLocs
%
% This script was created by Joe Cursons at the University of Melbourne
%   Systems Biology Laboratory:
%       joseph.cursons@unimelb.edu.au
%
% Last Updated: 24/09/15
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

%move through all sub-cellular localisations examined within these data
arrayLocStrings = { [ 'C' ];
                    [ 'N' ];
                    [ 'M' ] };
numLocalisations = length(arrayLocStrings);   

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
stringFormatForSampleLocationDataToLoad = 'MATLAB';
stringFormatForSampleLocationDataToWrite = 'TIFF';

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


 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Execute the required functions
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
 
%check that folder/path allocations are correct
% NB: I'm fairly sure exist(<>,'dir') doesn't work properly on UNIX - it
%  will always pass even when the folder is missing, I will update this
%  with an alternative function call when I find one
if (exist(strProcImgDataPath, 'dir') == 7) && (exist(strProcImgDataPath, 'dir') == 7),
    
    %move through all specified sub-cellular localisations
    for iLoc = 1:numLocalisations,
        %and move through all IF targets
        for iTarget = 1:numTargets,
            %and move through all patients
            for iPatient = 1:numPatients,
                
                %use these procedural loops and associated matrices to
                % specify relative folder paths
                numPatient = arrayPatients(iPatient);
                stringSegImgPath = [ strProcImgDataPath arrayTargetPaths{iTarget} strFoldSep 'Pat_' num2str(iPatient) strFoldSep ];
                stringSmpLocPath = [ strProcImgDataPath arrayTargetPaths{iTarget} strFoldSep 'Pat_' num2str(iPatient) strFoldSep arrayLocStrings{iLoc} strFoldSep ];
                
                %check that the specified sub-cellular localisation exists
                % for processing
                if (exist(stringSmpLocPath, 'dir') == 7),

                    %load the sample location data from the MATLAB format
                    % saved file
                    structSampleLocs = loadSampLocs(stringSmpLocPath, stringFormatForSampleLocationDataToLoad);
                    %then save the output as the specified format
                    % (binarised .TIFF files)
                    flagResult = saveSampLocs(structSampleLocs, stringSmpLocPath, stringFormatForSampleLocationDataToWrite);
                    if ~(flagResult == 1),
                        disp('warning: the extracted sample data do not appear to have been saved');
                    end
                    
                else
                    %throw back a warning - for a lot of the targets there
                    % are no plasma membrane data etc, so this will spam
                    % the command window
                    % To Do: specify output with different log levels 
                    disp(['warning: ' arrayTargetPaths{iTarget} ', Patient ' num2str(iPatient) ' has no data for localisation ' arrayLocStrings{iLoc}]);
                end

            end
        end
    end
    
else
    %throw back an error 
    disp(['error: data files cannot be located, please check the relative folder path for ' strProcImgDataPath ] );
    
end

