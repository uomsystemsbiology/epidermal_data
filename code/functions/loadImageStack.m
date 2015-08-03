function [ imagestack_multiframe ] = loadImageStack( path, stackname, n_images, channel )
%load_stack :: A function for reading an image stack into the memory
% [ imagestack_multiframe ] = load_stack( path, stackname, n_images, channel )
% ## Inputs:
% path - Full system path for the image stack, including the final "\"
% stackname - The name and extension format for the desired images. Include
% a * for the final z-index digit, and a # for the final channel-index
% digit.
% n_images - The number of images in the stack
% channel - The desired channel to load from the stack
% ## Outputs
% imagestack_multiframe - A multiframe array containing the image stack
%
% eg. im_stack = load_stack('F:\shared\full_share\confocal_results\08_08(aug)_13th\6_-_(Cx43+A555)\', '6_Series007_z00*_ch0#.tif', 47, 0);
%
% Created by Joe Cursons - j.cursons@auckland.ac.nz
% Last edited 24/03/09


% Identify the position of the * (identifies z-position marker) and #
% (identifies channel marker) characters within the stackname
for character_count = 1:length(stackname),
    if stackname(character_count) == '*',
        z_placer = character_count;
    elseif stackname(character_count) == '#',
        channel_placer = character_count;
    end
end

% Alter the stackname to include the z-position and channel
for image_count = 1:n_images,
    % NB: base 0 in file names, base 1 in arrays
    image_name = stackname;
    % Insert the image_number within the z-stack
    if length(int2str(image_count-1)) == 1,
        image_name(z_placer) = int2str(image_count-1);
    elseif length(int2str(image_count-1)) == 2,
        image_name(z_placer-1:z_placer) = int2str(image_count-1);
    elseif length(int2str(image_count-1)) == 3,
        image_name(z_placer-2:z_placer) = int2str(image_count-1);
    end
    % Insert the channel number
    image_name(channel_placer) = int2str(channel);
    imagestack(image_count).filename = image_name;
end

% Load all of the images according to the generated filenames and given
% path
for image_count = 1:n_images,
    imagestack(image_count).image = imread([path imagestack(image_count).filename]);
end

%% Convert all of the images into a multiframe image for display
imagestack_image_str = '';
for image_count = 1:n_images,
    imagestack_image_str = strcat(imagestack_image_str, 'imagestack(', int2str(image_count), ').image');
    if image_count < n_images,
        imagestack_image_str = strcat(imagestack_image_str, ', ');
    end
end

imagestack_multiframe = eval( [ 'cat(3,' imagestack_image_str ')' ] );