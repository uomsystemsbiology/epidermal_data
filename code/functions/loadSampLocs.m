function structOutSampleLocs = loadSampLocs( stringInDataPath, stringInDataFormat )
%% structOutSampleLocs = loadSampLocs( stringInDataPath, stringInDataFormat )
% This function takes a specified folder path for the sample location data,
%  and a string specifying the format, and attempts to load these data for
%  further analysis.
% 
%  Inputs:
%   - stringInDataPath: a string which contains the full path for the
%           specified data (i.e. including target/patient)
%   - stringInDataFormat: a string specifying the data format to load
%           
%  Output:
%   - structOutSampleLocs: a structured array contains sample locations for 
%           this patient/image data. Note that the SmpCent element has a
%           size (m*3), and the total number of samples for this
%           patient/image is m*length(structOutSampleLocs). m is usually 4
%           for older data sampled using scripts created by JC. The output
%           from this script when loading the image data has m = 1.
%       - SmpCent: an (m*3) array containing the (x,y,z) co-ordinates for 
%           m samples (as rows).
%       - ZPosition: an array specifying the Z-Positions sampled across the
%           full data set (used for the function which loads a sparse
%           arrays for the image data)
%
%  MATLAB Toolbox Dependencies:
%   - Image Processing Toolbox: this script directly calls imread
%                                    
%
% This MATLAB function has been released for the GigaScience Data Note:
%   Cursons et al. (2015). Spatially-transformed fluorescence image data 
%    for ERK-MAPK and selected proteins within human epidermis.
%    GigaScience. Submitted Sept 2015.
%   doi: not-yet-known
% 
% A more detailed description of the normalised distance co-ordinate can be
%  found in:
%   Cursons et al. (2015). Regulation of ERK-MAPK signaling in human 
%    epidermis. BMC Systems Biology. 
%   doi: 10.1186/s12918-015-0187-6
%
% This function was created by Joe Cursons:
%   joseph.cursons@unimelb.edu.au
%
% Last Updated: 24/09/15
%
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% User specified settings
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%specify an array of strings for legitimate data types
arrayAcceptedMATLABFormatStrings = {'MATLAB'; '.mat'; 'matlab'};
arrayAcceptedImageFormatStrings = {'TIFF'; '.TIFF'; 'tif'; '.tif'; 'tiff'; '.tiff'; 'image'; 'Image'};

%specify the string prefix for the sample location image files
stringSampleLocImagePrefix = 'im_sampleLocs_z';

 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Perform pre-processing
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  

%To do: write some checks that look for the specified data format and throw
%   back a warning/request user input on whether alternative data formats
%   are acceptable.
 
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Load the data and return
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
 
%check that the data format string is a legitimate form;
if any(strncmp(stringInDataFormat, arrayAcceptedMATLABFormatStrings, length(stringInDataFormat))),
    %load the MATLAB sample location files, which contains a 'SampleOutput'
    % struct
    load([ stringInDataPath 'sample_loc.mat']);
    
    %save to the specified output structure    
    structOutSampleLocs = SampleOutput;
    
elseif any(strncmp(stringInDataFormat, arrayAcceptedImageFormatStrings, length(stringInDataFormat))),
    %load the binarised TIFF files containing the sample location data, and
    % save this into an appropriate structured array
    
    %check the folder contents
    arrayFolderContents = dir(stringInDataPath);
    %determine which files correspond to sample location data
    arraySampLocFlag = false(length(arrayFolderContents), 1);
    for iFile = 1:length(arrayFolderContents),
        if strncmp(stringSampleLocImagePrefix, arrayFolderContents(iFile).name, length(stringSampleLocImagePrefix)),
            arraySampLocFlag(iFile) = true;
        end
    end
    %identify the file indices for subsequent processing
    arraySampLocIndices = find(arraySampLocFlag);
    numSampLocFiles = sum(arraySampLocFlag);
    
    %determine the total number of samples for memory allocation
    numTotalSamples = 0;
    for iZPos = 1:numSampLocFiles,
        %load the binarised TIFF file containing the sample locations
        stringSampleLocationFileName = arrayFolderContents(arraySampLocIndices(iZPos)).name;
        imZPosSamples = imread([stringInDataPath stringSampleLocationFileName]);
        %determine the number of samples from the length of the array for
        % non-zero entries (i.e. "1" within the binarised array, specifying
        % a samlpe location)
        numTotalSamples = numTotalSamples + length(find(imZPosSamples));
    end
    
    %create the output structured array
    structOutSampleLocs = struct('SmpCent', cell(numTotalSamples,1),  'ZPosition', cell(numTotalSamples,1));
    %initialise an output index counter
    numOutSample = 1;
    for iZPos = 1:numSampLocFiles,
        %load the binarised TIFF file containing the sample locations
        stringSampleLocationFileName = arrayFolderContents(arraySampLocIndices(iZPos)).name;
        imZPosSamples = imread([stringInDataPath stringSampleLocationFileName]);
        
        %and determine the "z" value
        numFileNameSuffixPosition = strfind(stringSampleLocationFileName, '.tiff');
        numFileNameZPosition = strfind(stringSampleLocationFileName, '_z');
        stringZPos = stringSampleLocationFileName((numFileNameZPosition+2):(numFileNameSuffixPosition-1));
        numZPos = uint8(str2double(stringZPos));
        
        %extract all of the sample locations
        [ arraySmpRow, arraySmpCol ] = find(imZPosSamples);
        
        %move through all samples and populate the output arrays
        for iSample = 1:length(arraySmpRow),
            
            %NB: the [row,col] output from find is reversed from the [x,y]
            % co-ordinates expected
            structOutSampleLocs(numOutSample).SmpCent = [ uint16(arraySmpCol(iSample)) uint16(arraySmpRow(iSample)) uint16(numZPos) ];
            structOutSampleLocs(numOutSample).ZPosition = numZPos;
            
            %increase the output index counter
            numOutSample = numOutSample + 1;
        end
        
    end
    
    
else
    disp('error: the specified data format could not be processed');
end
        
                
 
end

