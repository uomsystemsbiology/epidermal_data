function arrayOutLoessData = loadLoessCurve(stringInFolderPath, stringInDataFormat)
%% arrayOutLoessCurve = loadLoessCurve(stringInFolderPath, stringInDataFormat)
% This function takes 
% 
%  Inputs:
%   - stringInFolderPath: a string which contains the full path for the
%           folder containing the data to load (i.e. including target
%           /patient/localisation))
%   - stringInDataFormat: a string specifying the desired data format to
%           load
%
%  Output:
%   - arrayOutLoessCurve: a
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
%% Load the data and return
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  

    %check that the data format string is a legitimate form
    if any(strncmp(stringInDataFormat, arrayAcceptedMATLABFormatStrings, length(stringInDataFormat))),
        %check that the file exists
        if (exist([ stringInFolderPath 'loess_data.mat'], 'file') == 2),
            load([ stringInFolderPath 'loess_data.mat']);
            arrayOutLoessData = arrayLoessData;
        else
            disp(['error: the specified MATLAB loess data file could not be found in ' stringInFolderPath]);
            arrayOutLoessData = -1;
        end
        
    else
        disp('error: the specified data format could not be processed');
    end

end

