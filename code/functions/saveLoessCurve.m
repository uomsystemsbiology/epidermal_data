function flagResult = saveLoessCurve(arrayInLoessData, stringInFolderPath, stringInDataFormat)
%% flagResult = saveLoessCurve(arrayLoessCurve, stringFormatForLoessCurveDataToSave)
% This function takes 
% 
%  Inputs:
%   - arrayInLoessData: 
%   - stringInFolderPath: a string which contains the full path for the
%           folder to save the data (i.e. including target/patient
%           /localisation))
%   - stringInDataFormat: a string specifying the desired data format to
%           save
%
%  Output:
%   - flagResult: a
%
%  MATLAB Toolbox Dependencies:
%   - 
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
% Last Updated: 09/09/15
%
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% User specified settings
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%specify an array of strings for legitimate data types
arrayAcceptedMATLABFormatStrings = {'MATLAB'; '.mat'; 'matlab'};


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

%         %confirm that the user wishes to overwrite the data, as this can be
%         % destructive without an appropriate backup/version control system
%         strWarningInput = '';
%         while (~strncmp(strWarningInput, 'Y', 1) || ~strncmp(strWarningInput, 'y', 1) || ~strncmp(strWarningInput, 'N', 1) || ~strncmp(strWarningInput, 'n', 1))
%             strWarningInput = input(['Are you certain that you want to overwrite the SampleAnalysis file in ' stringInFolderPath '? [y/Y/n/N] '], 's');
%         end
% 
%         %if confirmed, save the data
%         if strncmp(strWarningInput, 'Y', 1) || strncmp(strWarningInput, 'y', 1),
            arrayLoessData = arrayInLoessData;
            save([ stringInFolderPath 'loess_data.mat'], 'arrayLoessData');
            flagResult = 1;
%         end
    else
        flagResult = 0;
        disp('error: the specified data format could not be processed');
    end

end

