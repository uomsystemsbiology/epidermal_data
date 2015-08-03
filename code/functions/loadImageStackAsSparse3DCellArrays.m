function [ arrayOutCellImageStack, arrayOutCellImageLayers, arrayOutCellImageBasalLamina, ...
    arrayOutCellImageBoundaryOne, arrayOutCellImageBoundaryTwo, ...
    arrayInCellImageOuterBoundary ] = loadImageStackAsSparse3DCellArrays( arrayInZSlices, ...
    stringInImageDirectory, stringInIndexImageFolder )
%load_3D_cell_arrays :: take the array of Z slices, and appropriate folder
%paths, and load the required cell arrays (empty on the Z-slices not
%required to save memory)
%
%Input:
%Output:

    % ======================= User Defined Settings =======================
    
    
    % ===================== Array/Variable Management =====================
    %variables based upon the input arrays
    numZSlicesSampled = size(arrayInZSlices,1);
    
    arrayDirContents = dir(stringInImageDirectory);
    numFilesInDir = size(arrayDirContents,1);
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
    
    %load the filenames into an array for searching        
    arrayFileNamesImageFolder = ls(stringInImageDirectory);
    numFilesinImageFolder = size(arrayFileNamesImageFolder,1);
    numImagesInStack = 0;
    for iFile = 1:numFilesinImageFolder,
        %search against the stack name up to the z-value as a stack-identifier
        if strfind(arrayFileNamesImageFolder(iFile,:), stringStackName(1:indexZSliceLoc+2));
            if strfind(arrayFileNamesImageFolder(iFile,:), ['ch0' num2str(numChannel)]),
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
    
    
    % ==================== Populate the Output Arrays =====================
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

