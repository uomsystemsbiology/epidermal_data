function [ arrayOutXCentre, arrayOutXBounds, arrayBoundaryIndices ] = calculateSpatialDivisions( numBins, arrayNormDistRatio )
%calculateSpatialDivisions
% a function that takes a specified number of spatial bins, and an array
% which describes the ratio of thickness for each tissue layer, then uses
% this information to calculate the edges of the spatial bins upon the
% normalised distance co-ordinate
% inputs:
%   - numBins: the total number of spatial bins (should be a multiple of
%               the summed arrayNormDistRatio)
%   - arrayNormDistRatio: an array containing the relative number of cells
%               within each tissue layer (to produce the total numBins)

    %create the required output arrays
    arrayBoundaryIndices = int16(zeros(4,1));
    arrayOutXCentre = double(zeros(numBins, 1));
    arrayOutXBounds = double(zeros(numBins+1, 1));

    %calculate any intermediate variables
    numNormDistFactor = sum(arrayNormDistRatio);
    numRatioBinsToFactor = numBins/numNormDistFactor;

    %identify the whole number (boundary) associated indices
    arrayBoundaryIndices(1) = 1;
    arrayBoundaryIndices(2) = (numRatioBinsToFactor*arrayNormDistRatio(1)) + 1;
    arrayBoundaryIndices(3) = numRatioBinsToFactor*(arrayNormDistRatio(1)+arrayNormDistRatio(2)) + 1;
    arrayBoundaryIndices(4) = numBins+1;

    arrayOutXBounds(  arrayBoundaryIndices(1):arrayBoundaryIndices(2)  ) = linspace(0, 1, (numRatioBinsToFactor*arrayNormDistRatio(1))+1 )';
    arrayOutXBounds(  arrayBoundaryIndices(2):arrayBoundaryIndices(3)  ) = linspace(1, 2, (numRatioBinsToFactor*arrayNormDistRatio(2))+1 )';
    arrayOutXBounds(  arrayBoundaryIndices(3):arrayBoundaryIndices(4)  ) = linspace(2, 3, (numRatioBinsToFactor*arrayNormDistRatio(3))+1 )';

    for iDivision = 1:numBins,
        arrayOutXCentre(iDivision) = ( arrayOutXBounds(iDivision) + arrayOutXBounds(iDivision+1) )/2;
    end

end

