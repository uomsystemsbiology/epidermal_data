function [ arrayLoessCurve ] = calculateLoessCurve( arrayXToSmooth, arrayYToSmooth, numSpatialBins, arrayNormDistRatio, numLoessWindowSize )
%calculateLoessCurve 
%   Detailed explanation goes here

    [ ~, ~, arrayBoundaryIndices ] = calculateSpatialDivisions( numSpatialBins, arrayNormDistRatio );

    numLayersInTissue = size(arrayBoundaryIndices,1)-1;
    arrayLayerBounds = arrayBoundaryIndices-1;

    arrayNormalisedXToSmooth = rescaleNormDist( arrayXToSmooth, numSpatialBins, arrayNormDistRatio );
    
    %sort the data for loess smoothing
    [arraySortedX, arraySortingPointer] = sort(arrayNormalisedXToSmooth);
    arrayTempYValuesSorted = arrayYToSmooth(arraySortingPointer);
    arraySmoothedY = cell(numLayersInTissue,1);
    disp('Performing loess smoothing, this may require extended processing times depending on the number of data points');
    for iTissueLayer = 1:numLayersInTissue,
        arraySamplesInTissueLayer = find(arraySortedX >= arrayLayerBounds(iTissueLayer) & ...
                                         arraySortedX <= arrayLayerBounds(iTissueLayer+1));
        arraySmoothedY{iTissueLayer} = smooth(arraySortedX(arraySamplesInTissueLayer),double(arrayTempYValuesSorted(arraySamplesInTissueLayer)),numLoessWindowSize,'rlowess');
    end
    arraySmoothedYCombined = cat(1, arraySmoothedY{1}, arraySmoothedY{2}, arraySmoothedY{3});

    %each point has a number of samples associated with it, but a unique loess
    %smoothed value (ie. there are heaps of repeats in the output arrays)
    numUniqueXPositions = length(unique(arraySortedX));


    %export the Loess curve data points in an 2xN array with the x- (row 1)
    % and y- (row 2) coordinates
    arrayLoessCurve = zeros(2,numUniqueXPositions, 'double');
    arrayLoessCurve(1,:) = unique(arraySortedX);
    for iXIndex = 1:numUniqueXPositions,
        %search for unique X values, and the associated Y value
        numSampleIndex = find(arraySortedX == arrayLoessCurve(1,iXIndex), 1);
        arrayLoessCurve(2,iXIndex) = arraySmoothedYCombined(numSampleIndex);
    end

    
    
end

