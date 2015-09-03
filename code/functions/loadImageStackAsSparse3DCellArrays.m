function [ arrayOutCellImageStack, arrayOutCellImageLayers, arrayOutCellImageBasalLamina, ...
    arrayOutCellImageBoundaryOne, arrayOutCellImageBoundaryTwo, ...
    arrayInCellImageOuterBoundary ] = loadImageStackAsSparse3DCellArrays( arrayInZSlices, ...
    stringInImageDirectory, stringInIndexImageFolder )
%% [ arrayOutCellImageStack, arrayOutCellImageLayers, arrayOutCellImageBasalLamina, ...
%    arrayOutCellImageBoundaryOne, arrayOutCellImageBoundaryTwo, ...
%    arrayInCellImageOuterBoundary ] 
%       = loadImageStackAsSparse3DCellArrays( arrayInZSlices, ...
%               stringInImageDirectory, stringInIndexImageFolder )
% This function is given folders containing image data and processed data,
%  and an array of desired z-positions, and the specified location is
%  extracted into sparse (non-sampled z-positions are empty) cell array.
%
%  Inputs:
%   - arrayInZSlices:
%   - stringInImageDirectory: full system path for the image stack, 
%           including the final path separator
%   - stringInIndexImageFolder: full system path for the segmentation data, 
%           including the final path separator
%
%  Output:
%   - arrayOutCellImageStack: a cell array containing the image data at
%           specified z-positions
%   - arrayOutCellImageLayers: a cell array containing the 
%           im_layers_z<#>.tif segmentation data  at specified z-positions
%   - arrayOutCellImageBasalLamina: a cell array containing the  
%           im_basal_z<#>.tiff segmentation data at specified z-positions
%   - arrayOutCellImageBoundaryOne: a cell array containing the  
%           im_bound1_z<#>.tiff segmentation data at specified z-positions
%   - arrayOutCellImageBoundaryTwo: a cell array containing the  
%           im_bound2_z<#>.tiff segmentation data at specified z-positions
%   - arrayInCellImageOuterBoundary:  a cell array containing the  
%           im_outer_bound_z<#>.tiff segmentation data at specified
%           z-positions
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
    %determine the total number of z-slices sampled
    numZSlicesSampled = size(arrayInZSlices,1);
    
    %extract the directory contents into a structured array for indexing
    arrayDirContents = dir(stringInImageDirectory);
    numFilesInDir = size(arrayDirContents,1);
    %move through the files and match to z-stack image data through the
    % 'Series' string
    iFile = 1;
    flagExitWhileLoop = 0;
    while (iFile <= numFilesInDir) && ~flagExitWhileLoop,
        if strfind(arrayDirContents(iFile).name, 'Series'),
            flagExitWhileLoop = 1;
        else
            iFile = iFile + 1;
        end
    end
    stringImageName = arrayDirContents(iFile).name;
    %search for the z-value identifier in the image name
    indexZSliceLoc = strfind(stringImageName, '_z');
    %determine the stack path as required by the load_stack function
    stringStackName = stringImageName;
    stringStackName(indexZSliceLoc+2:indexZSliceLoc+4) = '00*';
    indexChannelLoc = strfind(stringStackName, '_ch');
    stringStackName(indexChannelLoc+3:indexChannelLoc+4) = '0#';
    stringChannel = stringImageName(indexChannelLoc+3:indexChannelLoc+4);
    numChannel = uint16(str2double(stringChannel));
    
    numFilesinImageFolder = length(arrayDirContents);
    numImagesInStack = 0;
    for iFile = 1:numFilesinImageFolder,
        %search against the stack name up to the z-value as a stack-identifier
        if strfind(arrayDirContents(iFile).name, stringStackName(1:indexZSliceLoc+2));
            if strfind(arrayDirContents(iFile).name, ['ch0' num2str(numChannel)]),
                numImagesInStack = numImagesInStack + 1;
            end
        end

    end
        
    %load the image stack
    arrayImageStack = loadImageStack( stringInImageDirectory, stringStackName, numImagesInStack, numChannel );
    
    %create the output cell arrays
    arrayOutCellImageStack = cell(numImagesInStack,1);
    arrayOutCellImageLayers = cell(numImagesInStack,1);
    arrayOutCellImageBasalLamina = cell(numImagesInStack,1);
    arrayOutCellImageBoundaryOne = cell(numImagesInStack,1);
    arrayOutCellImageBoundaryTwo = cell(numImagesInStack,1);
    arrayInCellImageOuterBoundary = cell(numImagesInStack,1);
    
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Populate the Output Arrays
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  

    for iSampledSlice = 1:numZSlicesSampled,
        numZSlice = arrayInZSlices(iSampledSlice);
        arrayOutCellImageStack{numZSlice} = arrayImageStack(:,:,numZSlice);
        arrayOutCellImageLayers{numZSlice} = imread([stringInIndexImageFolder 'im_layers_z' num2str(numZSlice) '.tiff'],'tiff');
        arrayOutCellImageBasalLamina{numZSlice} = imread([stringInIndexImageFolder 'im_basal_memb_z' num2str(numZSlice) '.tiff'],'tiff');
        arrayOutCellImageBoundaryOne{numZSlice} = imread([stringInIndexImageFolder 'im_bound1_z' num2str(numZSlice) '.tiff'],'tiff');
        arrayOutCellImageBoundaryTwo{numZSlice} = imread([stringInIndexImageFolder 'im_bound2_z' num2str(numZSlice) '.tiff'],'tiff');
        arrayInCellImageOuterBoundary{numZSlice} = imread([stringInIndexImageFolder 'im_outer_bound_z' num2str(numZSlice) '.tiff'],'tiff');
    end
    

end

