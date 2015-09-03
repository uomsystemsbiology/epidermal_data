function [ numPixelsPerSample, arrayXOffset, arrayYOffset ] = initSamplingKernel( arraySamplingKernel )
%% [ numPixelsPerSample, arrayXOffset, arrayYOffset ] = initSamplingKernel( arraySamplingKernel )
% This function that takes an array specifying a sampling kernel (i.e. 
%  which pixels to sample around the central pixel)
%
%  Inputs:
%   - arraySamplingKernel: a binary array labelling the pixels which should
%           be sampled (around the central pixel for each specified region
%           of interest) for the signal intensity data
%
%  Output:
%   - numPixelsPerSample: the total number of pixels that will be sampled
%           from each region of interest
%   - arrayXOffset: a vector describing the x-offset to apply around the
%           central pixel when sampling
%   - arrayYOffset: a vector describing the y-offset to apply around the
%           central pixel when sampling
%
%  MATLAB Toolbox Dependencies:
%   - none
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
% Last Updated: 03/09/15
%
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Perform Pre-Processing
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  

    %determine the total number of pixels within the sampling kernel
    numPixelsPerSample = sum(double(arraySamplingKernel(:)));

    %determine the central pixel (hopefully within a kernel with an
    % odd-number of pixels along each edge)
    if iseven(size(arraySamplingKernel,1)),
        disp('warning: the function initSamplingKernel is designed to take a sampling kernel with an odd number of pixels in length (e.g. 3, 5, 7)');
    end
    numKernelCentre = ceil(size(arraySamplingKernel,1)/2);
 
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Produce the Output Arrays
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  

    %move through each pixel to be sampled and determine the x- and
    % y-offset
    arrayXOffset = zeros(numPixelsPerSample, 1, 'int16');
    arrayYOffset = zeros(numPixelsPerSample, 1, 'int16');
    [arrayRowPosition, arrayColumnPosition] = find(arraySamplingKernel);
    for iPixel = 1:numPixelsPerSample,
        arrayYOffset(iPixel) = arrayRowPosition(iPixel) - numKernelCentre;
        arrayXOffset(iPixel) = arrayColumnPosition(iPixel) - numKernelCentre;
    end

end

