function [ arrayMultiFrameImageStack ] = loadImageStack( stringFolderPath, stringStackName, numImages, numChannel )
%% arrayMultiFrameImageStack = loadImageStack( stringFolderPath, stringStackName, numImages, numChannel )
% This function reads an image stack into the memory as a large array for
%  subsequent processing.
%
%  Inputs:
%   - stringFolderPath: full system path for the image stack, including the
%           final path separator
%   - stringStackName: the name and extension format for the desired 
%           images. Include % a * for the final z-index digit, and a # for
%           the final channel-index digit.
%           e.g. '6_Series007_z00*_ch0#.tif'
%   - numImages: the number of images in the stack to be loaded
%   - numChannel: the number of the channel to load (usually 0 for the
%           single-target data)
%
%  Output:
%   - arrayMultiFrameImageStack: a multiframe (1024*1024*numImages) array
%           containing the image stack
%
%  MATLAB Toolbox Dependencies:
%   - none                      
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

    %identify the position of the * (z-position marker) and # (channel
    % marker) characters within the stackname
    for iChar = 1:length(stringStackName),
        if stringStackName(iChar) == '*',
            numZMarkerPos = iChar;
        elseif stringStackName(iChar) == '#',
            numChanMarkerPos = iChar;
        end
    end

    % Alter the stackname to include the z-position and channel
    for iImage = 1:numImages,
        % NB: base 0 in file names, base 1 in arrays
        stringImageName = stringStackName;
        %insert the image_number within the z-stack
        if length(int2str(iImage-1)) == 1,
            stringImageName(numZMarkerPos) = int2str(iImage-1);
        elseif length(int2str(iImage-1)) == 2,
            stringImageName(numZMarkerPos-1:numZMarkerPos) = int2str(iImage-1);
        elseif length(int2str(iImage-1)) == 3,
            stringImageName(numZMarkerPos-2:numZMarkerPos) = int2str(iImage-1);
        end
        %insert the channel number
        stringImageName(numChanMarkerPos) = int2str(numChannel);

        structImageStack(iImage).filename = stringImageName;
    end

    %load all images according to the generated filenames and given path
    for iImage = 1:numImages,
        structImageStack(iImage).image = imread([stringFolderPath structImageStack(iImage).filename]);
    end

%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Create the output array
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
    stringImageStackImage = '';
    for iImage = 1:numImages,
        stringImageStackImage = strcat(stringImageStackImage, 'structImageStack(', int2str(iImage), ').image');
        if iImage < numImages,
            stringImageStackImage = strcat(stringImageStackImage, ', ');
        end
    end

    arrayMultiFrameImageStack = eval( [ 'cat(3,' stringImageStackImage ')' ] );

end