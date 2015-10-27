function [ arrayImageOut ] = alterPixelIntensityForDisplay( arrayImageIn, numDarkLevel, numCoeff, numScalar )
% alter_pixel_intensity :: take an input image and apply a non-linear
%       transform of the pixel intensity:
%                   I_out = 255*numScalar*( (I_in/I_max)^gamma );
%
% inputs:
%   - arrayImageIn: a one-dimensional image array which needs to
%           undergo non-linear brightness modifications
%  	- numDarkLevel: the input image pixel intensity to be used for
%           thresholding background signal/noise
%               --> use n to bring all values <= n to 0
%               --> use 0 for no thresholding
%  	- numCoeff: the exponent used for the non-linear transform:
%               I_out = 255*numScalar*( (I_in/I_max)^gamma );
%           	--> use 0 to transform to a uniform distribution
%               --> use 0 < gamma < 1 to increase brightness
%               --> use 1 to perform no modification
%               --> use 1 < gamma to decrease brightness
%   - numScalar: the scalar value used to alter the maximum value/range of
%           pixel brightness
%               --> use < 1 to reduce total pixel brightness in the
%               specified channel
%               --> use == 1 for no change
%               --> use > 1 to increase total pixel brightness
%
% outputs:
%       - arrayImageOut: a one-dimensional image array which has undergone
%                   non-linear brightness modifications
%
% Created by Joe Cursons, The University of Melbourne, Australia
% email: joe ('dot') cursons ('at`) gmail (`dot`) com
% Last Modified 27/11/13


    
    %% ==================== Array/Variable Management =====================
    numMaxIntensity = double(max(arrayImageIn(:)));
    
    arrayImageOut = zeros(size(arrayImageIn,1), size(arrayImageIn,2), 'uint8');
    
    for YCoOrdinate = 1:size(arrayImageIn,1),
        for XCoOrdinate = 1:size(arrayImageIn,2),
            
            if (arrayImageIn(YCoOrdinate, XCoOrdinate) <= numDarkLevel),
                arrayImageOut(YCoOrdinate, XCoOrdinate) = uint8(0);
            else
                numPixelIntensity = double(arrayImageIn(YCoOrdinate, XCoOrdinate));
                arrayImageOut(YCoOrdinate, XCoOrdinate) = uint8(255*numScalar*((numPixelIntensity/numMaxIntensity)^numCoeff));
        end
    end

end

