function [ arrayOutXCentre, arrayOutXBounds, arrayBoundaryIndices ] = calculateSpatialDivisions( numBins, arrayNormDistRatio )
%% [ arrayOutXCentre, arrayOutXBounds, arrayBoundaryIndices ] = calculateSpatialDivisions( numBins, arrayNormDistRatio )
% This function is designed to take in a vector of relative tissue depths
%  (arrayNormDistRatio), together with a total number of spatial bins for
%  scaling; and then output the central positions (arrayOutXCentre) and
%  edges (arrayOutXBounds) of the resulting partitions, together with the
%  indices which specify the value at the layer boundaries, using the
%  normalised distance coordinate (for further details please refer to the
%  paper referenced below).
%
%  Inputs:
%   - numBins: the total number of spatial bins desired; note that this
%           should be a multiple of the summed arrayNormDistRatio
%   - arrayNormDistRatio: an array containing the relative number of cells
%           within each tissue layer; suggested value [1 4 2]
%
%  Output:
%   - arrayOutXCentre: an array containing the x-position for the centre of
%           each spatial partition
%   - arrayOutXBounds: an array containing the x-position for the edges of
%           each spatial partition
%   - arrayBoundaryIndices: an array containing the indices which
%           correspond to the layer boundaries
%   NB: it is recommended that you use the 'no output' value "~" for output
%       arrays which are not required
%
%  MATLAB Toolbox Dependencies:
%   - this function has no MATLAB toolbox dependencies
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
    %create the required output arrays
    arrayBoundaryIndices = int16(zeros(4,1));
    arrayOutXCentre = double(zeros(numBins, 1));
    arrayOutXBounds = double(zeros(numBins+1, 1));

    %calculate any intermediate variables
    numNormDistFactor = sum(arrayNormDistRatio);
    numRatioBinsToFactor = numBins/numNormDistFactor;
    
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Create the output arrays
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
 
    %identify the whole number (boundary) associated indices
    arrayBoundaryIndices(1) = 1;
    arrayBoundaryIndices(2) = (numRatioBinsToFactor*arrayNormDistRatio(1)) + 1;
    arrayBoundaryIndices(3) = numRatioBinsToFactor*(arrayNormDistRatio(1)+arrayNormDistRatio(2)) + 1;
    arrayBoundaryIndices(4) = numBins+1;

    %identify the x-position of boundaries for the partitions
    arrayOutXBounds(  arrayBoundaryIndices(1):arrayBoundaryIndices(2)  ) = linspace(0, 1, (numRatioBinsToFactor*arrayNormDistRatio(1))+1 )';
    arrayOutXBounds(  arrayBoundaryIndices(2):arrayBoundaryIndices(3)  ) = linspace(1, 2, (numRatioBinsToFactor*arrayNormDistRatio(2))+1 )';
    arrayOutXBounds(  arrayBoundaryIndices(3):arrayBoundaryIndices(4)  ) = linspace(2, 3, (numRatioBinsToFactor*arrayNormDistRatio(3))+1 )';

    %identify the x-position of centres for the partitions
    for iDivision = 1:numBins,
        arrayOutXCentre(iDivision) = ( arrayOutXBounds(iDivision) + arrayOutXBounds(iDivision+1) )/2;
    end

end

