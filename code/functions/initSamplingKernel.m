function [ numPixelsPerSample, arrayXOffset, arrayYOffset ] = initSamplingKernel( arraySamplingKernel )
%initSamplingKernel
% a function that takes an array specifying a sampling kernel (i.e. which
% pixels to sample aorund the central pixel)
% inputs:
%   - arraySamplingKernel: a binary array labelling the pixels which should
%           be sampled (around the central pixel for each specified region
%           of interest) for the signal intensity data
% outputs:
%   - numPixelsPerSample: the total number of pixels that will be sampled
%           from each region of interest
%   - arrayXOffset: a vector describing the x-offset to apply around the
%           central pixel when sampling
%   - arrayYOffset: a vector describing the y-offset to apply around the
%           central pixel when sampling


%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Produce the Output Arrays
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  

    numPixelsPerSample = sum(double(arraySamplingKernel(:)));

    %determine the central pixel (hopefully within a kernel with an
    %odd-number of pixels along each edge)
    numKernelCentre = ceil(size(arraySamplingKernel,1)/2);
    
    %move through each pixel to be sampled and determine the x- and
    %y-offset
    arrayXOffset = zeros(numPixelsPerSample, 1, 'int16');
    arrayYOffset = zeros(numPixelsPerSample, 1, 'int16');
    [arrayRowPosition, arrayColumnPosition] = find(arraySamplingKernel);
    for iPixel = 1:numPixelsPerSample,
        arrayYOffset(iPixel) = arrayRowPosition(iPixel) - numKernelCentre;
        arrayXOffset(iPixel) = arrayColumnPosition(iPixel) - numKernelCentre;
    end

end

