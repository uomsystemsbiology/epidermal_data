function [ arrayLoessCurve ] = calculateLoessCurve( arrayXToSmooth, arrayYToSmooth, numSpatialBins, arrayNormDistRatio, numLoessWindowSize )
%% arrayLoessCurve = calculateLoessCurve( arrayXToSmooth, arrayYToSmooth, numSpatialBins, arrayNormDistRatio, numLoessWindowSize )
% This function is designed to take in arrays containing data clouds
%  (y-data) against an x-coordinate (normalised distance) which is split 
%  into separate layers; see the papers referred below.
%
%  Inputs:
%   - arrayXToSmooth: a vector of x-coordinates, containing the location of
%           sampled pixels
%   - arrayYToSmooth: a vector of y-values at the specified x co-ordinates,
%           note that the script only expects a single observation at each
%           x position
%   - numSpatialBins: the number of bins for spatial partitioning of the
%           epidermis; note that this number should be a multiple for the
%           sum of arrayNormDistRatio (which by default sums to 7)
%   - arrayNormDistRatio: an array containing the relative number of cells
%           within each tissue layer; suggested value [1 4 2]
%   - numLoessWindowSize: a parameter for the loess smoothing (representing
%       the 'relative fraction' of the data cloud used for averaging;
%       recommended default: 0.5
%
%  Output:
%   - arrayLoessCurve: a (n,2) column vector with spatial positions
%   `   (normalised distance) in the first column, and loess-smoothed value
%       in the second column
%
%  MATLAB Toolbox Dependencies:
%   - Curve Fitting Toolbox: this function calls smooth
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

    %determine the indices which correspond to tissue layer boundaries,
    % dependent upon the specified number of spatial bins
    [ ~, ~, arrayBoundaryIndices ] = calculateSpatialDivisions( numSpatialBins, arrayNormDistRatio );
    numLayersInTissue = size(arrayBoundaryIndices,1)-1;
    arrayLayerBounds = arrayBoundaryIndices-1;

    %rescale the normalised distance value to the specified range
    arrayNormalisedXToSmooth = rescaleNormDist( arrayXToSmooth, numSpatialBins, arrayNormDistRatio );
    
    %sort the data prior to loess smoothing
    [arraySortedX, arraySortingPointer] = sort(arrayNormalisedXToSmooth);
    %apply the sorting index to the y-data
    arrayTempYValuesSorted = arrayYToSmooth(arraySortingPointer);
    %create the output cell array (each cell matches to a tissue layer)
    arraySmoothedY = cell(numLayersInTissue,1);
    
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Perform loess smoothing across the data
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
    disp('Performing loess smoothing, this may require extended processing times depending on the number of data points, please be patient..');
    
    %perform loess smoothing within each tissue layer individually to
    % minimise boundary effects
    for iTissueLayer = 1:numLayersInTissue,
        %identify samples within the specified tissue layer
        arraySamplesInTissueLayer = find(arraySortedX >= arrayLayerBounds(iTissueLayer) & ...
                                         arraySortedX <= arrayLayerBounds(iTissueLayer+1));
        %perform the loess smoothing over these samples
        arraySmoothedY{iTissueLayer} = smooth(arraySortedX(arraySamplesInTissueLayer),double(arrayTempYValuesSorted(arraySamplesInTissueLayer)),numLoessWindowSize,'rlowess');
    end
    %combine together the data across all tissue layers
    arraySmoothedYCombined = cat(1, arraySmoothedY{1}, arraySmoothedY{2}, arraySmoothedY{3});

    %each point has a number of samples associated with it, but a unique 
    % loesssmoothed value (ie. there are heaps of repeats in the output
    %  arrays)
    numUniqueXPositions = length(unique(arraySortedX));

 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Output in the Specified Format
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %      
    
    %export the loess curve data points in an 2xN array with the x- (row 1)
    % and y- (row 2) coordinates
    arrayLoessCurve = zeros(2,numUniqueXPositions, 'double');
    arrayLoessCurve(1,:) = unique(arraySortedX);
    for iXIndex = 1:numUniqueXPositions,
        %search for unique X values, and the associated Y value
        numSampleIndex = find(arraySortedX == arrayLoessCurve(1,iXIndex), 1);
        arrayLoessCurve(2,iXIndex) = arraySmoothedYCombined(numSampleIndex);
    end

end

