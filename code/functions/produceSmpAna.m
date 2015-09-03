function [ structSigDataOut ] = produceSmpAna( structSampleLocs, arraySamplingKernel, stringImgDataPath, stringSegImgPath )
%% structSigDataOut = produceSmpAna( structSampleLocs, arraySamplingKernel, stringImgDataPath, stringSegImgPath )
% This function which takes a structured array of sample
%  locations, a specified sampling kernel (for selecting pixels around the
%  region of interest) and the folder paths for image and segmentation 
%  data, and produces a structured array containing signal intensity and 
%  pixel location data for subsequent analysis
%
%  Inputs:
%   - structSampleLocs: a structured array containing 3D pixel co-ordinates
%           (x, y, z) within a 'SmpCent' Nx3 array with N sample locations
%   - arraySamplingKernel: an binarised array showing pixels to extract the
%           signal intensity data around the central pixel, for each sample
%           location
%           NB: that around line 91 a call to the unique function is made
%            to ensure that there are no over-lapping pixels being sampled
%            after applying the sampling kernel at the specified locations
%   - stringImgDataPath: the folder path for the image data stack
%   - stringSegImgPath: the folder path for the tissue segmentation and 
%           sample location data
%
%  Output:
%   - structSigDataOut: a structured array containing:
%   	- CoOrds: a vector with the x, y and z co-ordinates of the sampled
%       	pixel
%   	- Layer: an integer with the tissue layer
%   	- BasDist: the distance (in pixels) from the basement membrane
%   	- NormDist: layer-normalized distance (linear interpolation within 
%       	each tissue layer)
%   	- SigInt: the signal intensity (uint8) associated with the pixel
%
%  MATLAB Toolbox Dependencies:
%   - none (?)
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
%% Perform Pre-Processing
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  

    %depending upon the OS, use a different folder separater (forward vs
    % back slash)
    if ispc,
        strFoldSep = '\';
    elseif isunix,
        strFoldSep = '/';
    else
        disp(['warning: cannot determine the operating system, defaulting to ' ...
                'forward slash for the folder path separator' ] );
        strFoldSep = '/';
    end

    %determine properties of the sampling kernel
    [ numPixelsPerSample, arrayXOffset, arrayYOffset ] = initSamplingKernel( arraySamplingKernel );

    %extract information about the number of samples (ROIs) which will be
    %used for signal intensity data (with the associated sampling kernel)
    numObjectsSampled = length(structSampleLocs);
    numSamplesPerObject = size(structSampleLocs(1).SmpCent,1);
    
    %extract the co-ordinates of all sampled pixels
    arrayAllPixelCoOrds = zeros(numObjectsSampled*numSamplesPerObject*numPixelsPerSample,3,'uint16');
    for iObject = 1:numObjectsSampled,
        for iSample = 1:numSamplesPerObject,
            %determine the indices for all pixels associated with this object/ROI
            numStartIndex = (iObject-1)*numSamplesPerObject*numPixelsPerSample + (iSample-1)*numPixelsPerSample;
            for iPix = 1:numPixelsPerSample,
                %calculate the x-, y- and z-positions
                arrayAllPixelCoOrds(numStartIndex+iPix,1) = uint16(int16(structSampleLocs(iObject).SmpCent(iSample,1)) + int16(arrayXOffset(iPix)));
                arrayAllPixelCoOrds(numStartIndex+iPix,2) = uint16(int16(structSampleLocs(iObject).SmpCent(iSample,2)) + int16(arrayYOffset(iPix)));
                arrayAllPixelCoOrds(numStartIndex+iPix,3) = uint16(structSampleLocs(iObject).SmpCent(iSample,3)); 
            end
        end
    end

    %ensure the pixel co-ordinates are unique to prevent re-sampling
    arrayUniquePixelCoOrds = unique(arrayAllPixelCoOrds, 'rows');
    numUniquePixels = length(arrayUniquePixelCoOrds);


    %load the images required for sampling
    arrayUniqueZPositions = unique(arrayUniquePixelCoOrds(:,3));
    [ arrayImageData, arrayImageLayers, arrayImageBasalLaminas, arrayImageBoundaryOnes, arrayImageBoundaryTwos, arrayImageOuterBoundaries ] = loadImageStackAsSparse3DCellArrays( arrayUniqueZPositions, stringImgDataPath, stringSegImgPath );
    

    %for error output, determine the target protein
    arraySegImgPathSeps = strfind(stringSegImgPath, strFoldSep);
    stringTarget = stringSegImgPath((arraySegImgPathSeps(end-2)+1):(arraySegImgPathSeps(end-1)-1));

%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Populate the Output Arrays
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %      

    %initialise the output data structure
    structSigDataOut = struct('CoOrds', cell(numUniquePixels,1), 'Layer', cell(numUniquePixels,1), 'BasDist', cell(numUniquePixels,1), 'NormDist', cell(numUniquePixels,1), 'SigInt', cell(numUniquePixels,1));

    %move through each pixel
    for iPix = 1:numUniquePixels,
        %initialise vectors within the output data structure
        structSigDataOut(iPix).CoOrds = zeros(1,3,'uint16');
        structSigDataOut(iPix).Layer = zeros(1,1,'int16');
        structSigDataOut(iPix).BasDist = zeros(1,1,'double');
        structSigDataOut(iPix).NormDist = zeros(1,1,'double');
        structSigDataOut(iPix).SigInt = zeros(1,1,'uint8');
        
        %extract the pixel co-ordinates
        structSigDataOut(iPix).CoOrds = arrayAllPixelCoOrds(iPix,:);
        numX = structSigDataOut(iPix).CoOrds(1);
        numY = structSigDataOut(iPix).CoOrds(2);
        numZ = structSigDataOut(iPix).CoOrds(3);

        %read the layer information from the associated image stack
        structSigDataOut(iPix).Layer = int16(arrayImageLayers{numZ}(numY, numX));
        %read the signal intensity data from the associated image stack
        structSigDataOut(iPix).SigInt = uint8(arrayImageData{numZ}(numY,numX));

        %note that it is much quicker to pass large vectors of co-ordinates
        %into calculateDistancesToBoundary, so BasDist and NormDist are
        %calculated within an alternative loop (below)
    end
    
    %move through each sampled z-slice and calculate the normalised 
    %distance values against the associated layer/boundary images - note
    %that it is MUCH quicker to pass large arrays of numbers into
    %calculateDistancesToBoundary than to calculate each pixel
    %co-ordinate individually
    for iZPos = 1:length(arrayUniqueZPositions),
        numZPos = arrayUniqueZPositions(iZPos);
        arrayInLayerPointer = find(arrayUniquePixelCoOrds(:,3) == numZPos);

        array2DCoOrds = zeros(2,length(arrayInLayerPointer), 'uint16');
        for iPix = 1:length(arrayInLayerPointer),
            numPointer = arrayInLayerPointer(iPix);
            array2DCoOrds(:,iPix) = structSigDataOut(numPointer).CoOrds(1:2);
        end

        %pass to the calculateDistancesToBoundary function
        arrayBasalLaminaDistance = calculateDistancesToBoundary(arrayImageBasalLaminas{numZPos}, array2DCoOrds);
        arrayBoundaryOneDistance = calculateDistancesToBoundary(arrayImageBoundaryOnes{numZPos}, array2DCoOrds);
        arrayBoundaryTwoDistance = calculateDistancesToBoundary(arrayImageBoundaryTwos{numZPos}, array2DCoOrds);
        arrayOuterBoundaryDistance = calculateDistancesToBoundary(arrayImageOuterBoundaries{numZPos}, array2DCoOrds);

        %read back into the nrmalised distance and basal distance values
        for iPix = 1:length(arrayInLayerPointer),
            numPointer = arrayInLayerPointer(iPix);
            numBasDist = arrayBasalLaminaDistance(iPix);
            numBoundOneDist = arrayBoundaryOneDistance(iPix);
            numBoundTwoDist = arrayBoundaryTwoDistance(iPix);
            numOuterBoundDist = arrayOuterBoundaryDistance(iPix);
            
            structSigDataOut(numPointer).BasDist = numBasDist;

            if (structSigDataOut(numPointer).Layer == 1),
                structSigDataOut(numPointer).NormDist = 0 + numBasDist/(numBasDist + numBoundOneDist);
            elseif (structSigDataOut(numPointer).Layer == 2),
                structSigDataOut(numPointer).NormDist = 1 + numBoundOneDist/(numBoundOneDist + numBoundTwoDist);
            elseif (structSigDataOut(numPointer).Layer == 3),
                structSigDataOut(numPointer).NormDist = 2 + numBoundTwoDist/(numBoundTwoDist + numOuterBoundDist);           
            elseif (structSigDataOut(numPointer).Layer == -1),
                structSigDataOut(numPointer).NormDist = 0 - numBasDist/(numBasDist + numBoundOneDist);
            elseif (structSigDataOut(numPointer).Layer == -3),
                structSigDataOut(numPointer).NormDist = 3 + numOuterBoundDist/(numBoundTwoDist + numOuterBoundDist);    
            else
                numX = structSigDataOut(numPointer).CoOrds(1);
                numY = structSigDataOut(numPointer).CoOrds(2);
                numZ = structSigDataOut(numPointer).CoOrds(3);
                disp(['error: invalid value for the sample layer (x,y,z): (' num2str(numX) ', ' num2str(numY) ', ' num2str(numZ) ')']);
                disp(['       target: ' stringTarget ]);
                structSigDataOut(numPointer).NormDist = NaN;  
            end
        end
    end

end

