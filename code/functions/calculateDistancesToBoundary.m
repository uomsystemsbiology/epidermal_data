function [ arrayEuclidianDist ] = calculateDistancesToBoundary( imageEdge, arrayPixelCoOrds )
%% arrayEuclidianDist = calculateCellDistances(imageEdge, arrayPixelCoOrds)
% This function calculates the minimum distance for an array of pixel 
%   co-ordinates, to a continuous boundary.
%
%  Inputs:
%   - imageEdge: a binarised image file containing the tissue-layer 
%           boundary, from which the distance of specified x,y co-ordinates
%           are to be calculated
%   - arrayPixelCoOrds: a row/column vector containing x,y co-ordinates of
%           the positions for distance calculation
%  Output:
%   - arrayEuclidianDist: an array (double precision) containing the
%           Euclidian distance between the specified centroid (x,y) 
%           positions and their nearest pixel on the input imageEdge
%  MATLAB Toolbox dependencies:
%   - Image Processing Toolbox: this function calls bwlabel and regionprops
%           
%
% This MATLAB function has been released for the GigaScience Data Note:
%   Cursons et al. (2015). Spatially-transformed fluorescence image data 
%    for ERK-MAPK and selected proteins within human epidermis.
%    GigaScience. Submitted Sept 2015.
%    doi: not-yet-known
%
% This function was created by Joe Cursons:
%   joseph.cursons@unimelb.edu.au
%
% Last Updated: 03/09/15
%
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Check Input Format
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
    %if the binarised imageEdge has been created within image processing
    % software it may actually be a '3D' (RGB) .TIFF file; if so, convert
    % to a one-dimensional binarised array
    if size(imageEdge,3) > 1,
        imageEdgeBin = bwlabel( imageEdge(:,:,1) );
    else
        imageEdgeBin = bwlabel( imageEdge );
    end
    
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Perform Pre-Processing
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
    %use the regionprops function to extract pixel co-ordinates from the
    % binarised edge file
    structEdgeProps = regionprops( imageEdgeBin, 'PixelList' );
    arrayEdgePixelCoOrds = [ structEdgeProps(1).PixelList' ];
    numPixelEdgeCoOrds = size( arrayEdgePixelCoOrds, 2 );
    
    %transpose the PixelList - the regionprops function extracts a column
    % vector, the dist function requires a row vector
    if (size(arrayPixelCoOrds,2) == 2) && (size(arrayPixelCoOrds,1) > 1),
        arrayPixelCoOrds = arrayPixelCoOrds';
    end

    %append the centroid co-ordinates to the end of the edge co-ordinates
    arrayEdgePixelCoOrds = [ arrayEdgePixelCoOrds, arrayPixelCoOrds ];
    numTotalPix = size( arrayEdgePixelCoOrds, 2 );
    
    %use the dist function to calculate the distance between all pixels
    % within the resulting array
    arrayEdgeDistances = dist( double(arrayEdgePixelCoOrds) );

 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Output in the Specified Format
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
 
    %identify the shortest distance between the edge pixels and the
    % specified objects
    arrayEuclidianDist = min(arrayEdgeDistances(1:numPixelEdgeCoOrds, numPixelEdgeCoOrds+1:numTotalPix));

end

