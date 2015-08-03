function [ arrayMinEuclidDist ] = calculateDistancesToBoundary( imEdge, arrayPixelCoOrds )
%calculateDistancesToBoundary :: calculate the minimum distance for an
%array of pixel co-ordinates, to a continuous boundary
%   
%   This script moves takes an array of co-ordinates for specified pixels
%(can be a single value, but it is much quicker to submit a large array
%with multiple pixel co-ordinates than to loop and submit the co-ordinates 
%of individual pixels for calculation), and a binarised image containing a
%boundary, and calculates the shortest distance to the boundary for each
%specified pixel co-ordinate.
% Inputs:
%   - imEdge: a binarised image containing a boundary (e.g. a tissue layer
%           boundary)
%   - arrayPixelCoOrds: a 2xN array of pixel co-ordinates (x in row 1, y in
%           row 2) to determine the minimum distance
% Outputs:
%   - arrayMinEuclidDist: the minimum Euclidian distance (in pixels) from
%           each specified pixel to the specified boundary

    %ensure the boundary edge is binarised
    arrayEdgeLabelled = logical( imEdge );
    %and use the regionprops function to extract the pixel co-ordinates of
    % the boundary
    arrayEdgeProps = regionprops( arrayEdgeLabelled, 'PixelList' );

    %transpose the PixelList as regionprops function extracts a column
    % vector with each pixel corresponding to a row; while the dist 
    % function requires each pixel to correspond to a row
    arrayEdgePixCoOrds = arrayEdgeProps(1).PixelList';
    numEdgePix = size( arrayEdgePixCoOrds, 2 );

    %concatenate the imEdge pixel list with the specified pixel
    % co-ordinates
    arrayEdgePixCoOrds = [ arrayEdgePixCoOrds, arrayPixelCoOrds ];
    numTotalPix = size( arrayEdgePixCoOrds, 2 );
    %and call the intrinsic MATLAB 'dist' function
    arrayEdgeDist = dist( double(arrayEdgePixCoOrds) );

    %identify the minimum distance for the specified pixel co-ordinates (in
    % columns numEdgePix+1:numTotalPix)
    arrayMinEuclidDist = min(arrayEdgeDist(1:numEdgePix, numEdgePix+1:numTotalPix));

end

