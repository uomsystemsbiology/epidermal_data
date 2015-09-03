function [ arrayOutNormDistRescaled ] = rescaleNormDist( arrayInNormDist, numInSpatialBins, arrayInNormDistRatio )
%% arrayOutNormDistRescaled = rescaleNormDist( arrayInNormDist, numInSpatialBins, arrayInNormDistRatio )
% This function is designed to take an input array of normalised distance
%  values and rescale them for plotting according to to the number of
%  spatial bins, or general manipalation of normalised distance values
%  (rescaling the number of spatial bins etc)
% 
%  Inputs:
%   - arrayInNormDist: a vector of normalised distance values to be 
%           rescaled over the specified number of spatial bins
%   - numSpatialBins: the number of bins for spatial partitioning of the
%           epidermis; note that this number should be a multiple for the
%           sum of arrayNormDistRatio (which by default sums to 7)
%   - arrayInNormDistRatio: an array containing the relative number of 
%           cells within each tissue layer; suggested value [1 4 2]
%
%  Output:
%   - arrayOutNormDistRescaled: a
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
% Last Updated: 03/09/15
%
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Perform pre-processing
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
 
    %determine the total number of samples to be rescaled
    numSamples = size(arrayInNormDist,1);
    %and the relative number of spatial partitions per 'cell' from the
    % normalised distance ratio
    numSamplesPerRatio = numInSpatialBins/sum(arrayInNormDistRatio);
    
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Output the Specified Array
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %      
    %create the output array
    arrayOutNormDistRescaled = double(zeros(numSamples,1));
    
    %step through each sample
    for iSample = 1:numSamples,
        %identify the layer from the whole-number component
        numSampleNormDist = arrayInNormDist(iSample);
        numSampleLayer = floor(numSampleNormDist);
        
        %the fractional distance is then scaled according to the layer
        %distance
        numSampleFractionalDistance = numSampleNormDist-numSampleLayer;
        if numSampleLayer == 0,
            %separate the 'fraction' to be rescaled
            numSampleNormDistFractionRescaled = numSampleFractionalDistance*numSamplesPerRatio*arrayInNormDistRatio(1);
            %from the 'tissue layer integer' to be rescaled
            numSampleNormDistLayerRescaled = 0;
        elseif numSampleLayer == 1,
            %separate the 'fraction' to be rescaled
            numSampleNormDistFractionRescaled = numSampleFractionalDistance*numSamplesPerRatio*arrayInNormDistRatio(2);
            %from the 'tissue layer integer' to be rescaled
            numSampleNormDistLayerRescaled = numSamplesPerRatio*arrayInNormDistRatio(1);
        elseif numSampleLayer == 2,
            %separate the 'fraction' to be rescaled
            numSampleNormDistFractionRescaled = numSampleFractionalDistance*numSamplesPerRatio*arrayInNormDistRatio(3);
            %from the 'tissue layer integer' to be rescaled
            numSampleNormDistLayerRescaled = numSamplesPerRatio*sum(arrayInNormDistRatio(1:2));
        elseif numSampleLayer == 3,
            %separate the 'fraction' to be rescaled
            numSampleNormDistFractionRescaled = 0;
            %from the 'tissue layer integer' to be rescaled
            numSampleNormDistLayerRescaled = numSamplesPerRatio*sum(arrayInNormDistRatio(1:3));
        end
        
        %recombine the 'tissue layer integer' and 'fraction'
        numSampleNormDistRescaled = numSampleNormDistLayerRescaled + numSampleNormDistFractionRescaled;
        
        %output in the appropriate vector position
        arrayOutNormDistRescaled(iSample) = numSampleNormDistRescaled;
    end

end

