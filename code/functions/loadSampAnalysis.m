function structOutSampleAnalysis = loadSampAnalysis( stringInDataPath, stringInDataFormat )
%% structOutSampleAnalysis = loadSampAnalysis( stringInDataPath, stringInDataFormat )
% This function loads sample analysis data, using strings to specify the
%  folder path and format of the data to be loaded.
% 
%  Inputs:
%   - stringInDataPath: a string which contains the full path for the
%           specified data (i.e. including target/patient)
%   - stringInDataFormat: a string specifying the desired data format to
%           load
%
%  Output:
%   - structOutSampleAnalysis: a structured array which contains sample 
%           analytics for further analysis, with a length corresponding to
%           the total number of unique sampled pixels
%       - CoOrds: a vector of (x,y,z) co-ordinates from the image data
%       - Layer: an integer specifying the tissue layer of the data point
%       - BasDist: a double specifying the Euclidian distance (pixels) to
%           the nearest point of the basal layer
%       - NormDist: a double specifying the normalised distance value for
%           the data point
%       - SigInt: an unsigned 8-bit integer specifying the signal intensity
%           at the data point
%
%  MATLAB Toolbox Dependencies:
%   - this function has no dependencies upon MATLAB toolboxes
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

    %check that the data format string is a legitimate form;
    if any(strncmp(stringInDataFormat, arrayAcceptedMATLABFormatStrings, length(stringInDataFormat))),
        %load the MATLAB sample location files, which contains a
        % 'SampleAnalysis' struct
        load([ stringInDataPath 'sample_ana.mat']);
        %save to the specified output structure    
        structOutSampleAnalysis = SampleAnalysis;
    elseif any(strncmp(stringInDataFormat, arrayAcceptedMATLABFormatStrings, length(stringInDataFormat))),

    else
        disp('error: the specified data format could not be processed');
    end
        

end

