function flagResult = saveSampLocs( structInLocsToSave, stringInFolderPath, stringInDataFormat )
%% flagResult = saveSampAnalysis( stringInFolderPath, stringInDataFormat )
% This function takes a structured array containing sample analysis data,
%  together with strings specifying the path and format for saving the
%  data.
% 
%  Inputs:
%   - structInDataToSave: the structured array containing the sample
%           analysis data to be saved
%   - stringInFolderPath: a string which contains the full path for the
%           folder to save the data (i.e. including target/patient
%           /localisation)
%   - stringInDataFormat: a string specifying the desired data format to
%           save
%
%  Output:
%   - flagResult: a numerical value specifying whether the operation was
%           successfully performed or not (this may be modified in future
%           values to integrate with error logging)
%
%  MATLAB Toolbox Dependencies:
%   - Image Processing Toolbox: this script directly calls imwrite
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

 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Perform pre-processing
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  

%To do: write some checks that look for the specified data format and throw
%   back a warning/request user input on whether alternative data formats
%   are acceptable.

 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Save the data and return
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
 
%check that the data format string is a legitimate form;
if any(strncmp(stringInDataFormat, arrayAcceptedMATLABFormatStrings, length(stringInDataFormat))),
    %save the structured array containing the sample analysis to a MATLAB 
    % data file, as a structured array called 'SampleOutput'
    SampleOutput = structInLocsToSave;
    
    %confirm that the user wishes to overwrite the data, as this can be
    % destructive without an appropriate backup/version control system
    strWarningInput = 'X';
    while ~( strncmp(strWarningInput, 'Y', 1) || strncmp(strWarningInput, 'y', 1) || ...
            strncmp(strWarningInput, 'N', 1) || strncmp(strWarningInput, 'n', 1) ),
        strWarningInput = input(['Are you certain that you want to overwrite the SampleAnalysis file in ' stringInFolderPath '? [y/Y/n/N]'], 's');
    end
    
    %if confirmed, save the data
    if strncmp(strWarningInput, 'Y', 1) || strncmp(strWarningInput, 'y', 1),
        save([ stringInFolderPath 'sample_loc.mat'], 'SampleOutput');
        flagResult = 1;
    else
        disp(['warning: the data were not saved in the ' stringInDataFormat ' format, as permission was not given to overwrite the existing file']);
    end
    
elseif any(strncmp(stringInDataFormat, arrayAcceptedImageFormatStrings, length(stringInDataFormat))),
    %extract pixel locations and save to a binarised TIFF file 
    
    %determine how many sample locations there are for these IF data
    numObjectsSampled = length(structInLocsToSave);
    numSamplesPerObject = size(structInLocsToSave(1).SmpCent,1);
    
    %extract the sample location data - do the z-positions as a separate
    % array for subsequent indexing
    arraySampleZPos = zeros(numObjectsSampled*numSamplesPerObject, 1, 'uint8');
    arraySampleXYPos = zeros(numObjectsSampled*numSamplesPerObject, 2, 'uint16');
    for iObject = 1:numObjectsSampled,
        numOutputBase = (iObject-1)*numSamplesPerObject;
        arraySampleZPos((numOutputBase+1):(numOutputBase+numSamplesPerObject)) = uint8(structInLocsToSave(iObject).ZPosition);
        for iSample = 1:numSamplesPerObject,
            arraySampleXYPos((numOutputBase+iSample),:) = uint16(structInLocsToSave(iObject).SmpCent(iSample,1:2));
        end
    end
    
    %determine the unique z-positions so that corresponding sample location
    % TIFF files can be created for each one
    arrayUniqueZPositions = unique(arraySampleZPos);
    for iZPos = 1:length(arrayUniqueZPositions),
        %determine which samples correspond to this z-position
        numZPos = arrayUniqueZPositions(iZPos);
        arraySamplesInZPos = find(arraySampleZPos == numZPos);
        
        %create a binarised array
        imSampLocs = false(1024,1024);
        %populate it with the sample locations
        for iSample = 1:length(arraySamplesInZPos),
            numYPos = arraySampleXYPos(arraySamplesInZPos(iSample),2);
            numXPos = arraySampleXYPos(arraySamplesInZPos(iSample),1);
            imSampLocs(numYPos, numXPos) = true;
        end
                
        %confirm that the user wishes to overwrite the data, as this can be
        % destructive without an appropriate backup/version control system
        strWarningInput = 'X';
        while ~( strncmp(strWarningInput, 'Y', 1) || strncmp(strWarningInput, 'y', 1) || ...
                 strncmp(strWarningInput, 'N', 1) || strncmp(strWarningInput, 'n', 1) ),
            strWarningInput = input(['Are you certain that you want to overwrite the sample location TIFF files in ' stringInFolderPath ' [y/Y/n/N]: '], 's');
        end

        %if confirmed, save the data as a binarised tiff file
        if strncmp(strWarningInput, 'Y', 1) || strncmp(strWarningInput, 'y', 1),
            imwrite(imSampLocs, [ stringInFolderPath 'im_sampleLocs_z' num2str(numZPos,'%i') '.tiff'], 'TIFF');
            flagResult = 1;
        else
            disp(['warning: the data were not saved in the ' stringInDataFormat ' format, as permission was not given to overwrite the existing file']);
        end
        
    end
    
else
    disp('error: the specified data format could not be processed');
    flagResult = 0;
end
 
 

end

