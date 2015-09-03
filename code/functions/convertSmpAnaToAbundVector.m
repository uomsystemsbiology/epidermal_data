function [ arrAbundVect ] = convertSmpAnaToAbundVector( structSampleInfo, numSpatialBins, arrayNormDistRatio, numLoessWindowSize, flagPerformZNorm )
%% arrAbundVect = convertSmpAnaToAbundVector( structSampleInfo, numSpatialBins, arrayNormDistRatio, numLoessWindowSize, flagPerformZNorm )
% This function is designed to take an input data structure which contains
%  information on the sample data (structSampleInfo), together with spatial
%  partitioning information (numSpatialBins, arrayNormDistRatio), smoothing
%  settings (numLoessWindowSize) and output normalisation settings
%  (flagPerformZNorm)
%
%  Inputs:
%   - structSampleInfo: a structured array containing a information on
%           sampled pixels, such as the normalised distance position of the
%           sample, and corresponding pixel intensities
%   - numSpatialBins: the number of bins for spatial partitioning of the
%           epidermis; note that this number should be a multiple for the
%           sum of arrayNormDistRatio (which by default sums to 7)
%   - arrayNormDistRatio: an array containing the relative number of cells
%           within each tissue layer; suggested value [1 4 2]
%   - numLoessWindowSize: a parameter for the loess smoothing (representing
%           the 'relative fraction' of the data cloud used for averaging;
%           recommended default: 0.5
%   - flagPerformZNorm: a binarised flag as to whether the data should
%           undergo z-score normalisation (i.e. mapped to standard 
%           deviations around the mean)
%
%  Output:
%   - arrAbundVect: a vector of interpolated, smoothed sample intensity 
%           (potentially after z-score normalisation) at the centre of
%           specified spatial bins
%
%  MATLAB Toolbox Dependencies:
%   - Curve Fitting Toolbox: this function calls calculateLoessCurve, which
%           uses smooth
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

    %call the calculate spatial divisions function for the specified
    % normalised distance parameters
    [ arrayOutXCentres, ~, ~ ] = calculateSpatialDivisions( numSpatialBins, arrayNormDistRatio );
    
    %rescale the x co-ordinates accordingly to identify the required output
    % data points (at the centre of each spatial partition)
    arrayOutXCentresRescaled = rescaleNormDist( arrayOutXCentres, numSpatialBins, arrayNormDistRatio );

%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% populate required arrays
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  

    %examine the size of the input array
    numSamples = length(structSampleInfo);
    numPixelsPerSample = length(structSampleInfo(1).SigInt);

    %move the data into temporary arrays for further processing
    arrayTempSamplePositions = zeros(numSamples, 1, 'double');
    arrayTempSampleValues = zeros(numSamples, numPixelsPerSample, 'uint8');
    for iSample = 1:numSamples,
        arrayTempSamplePositions(iSample) = structSampleInfo(iSample).NormDist;
        arrayTempSampleValues(iSample, 1:numPixelsPerSample) = structSampleInfo(iSample).SigInt;
    end

    %if specified, perform z-score normalisation
    if flagPerformZNorm,
        %calculate the mean and standard deviation
        arraySampleMean = mean(double(arrayTempSampleValues(:)));
        arraySampleStDev = std(double(arrayTempSampleValues(:)));
        %centre around the mean and normalised against the standard
        % deviation
        arrayDataToSmooth = (double(arrayTempSampleValues) - arraySampleMean)/arraySampleStDev;
    else
        arrayDataToSmooth = arrayTempSampleValues;
    end
    
    %rescale the x co-ordinates of the sample data
    arrayRescaledSamplePositions = rescaleNormDist( arrayTempSamplePositions, numSpatialBins, arrayNormDistRatio );

    %perform loess smoothing
    arrayLoessCurveXY = calculateLoessCurve( arrayRescaledSamplePositions, arrayDataToSmooth, numSpatialBins, arrayNormDistRatio, numLoessWindowSize );
    
    %apply cubic hermite interpolation to extract the loess curve value at
    % the specified output positions
    arrAbundVect = interp1(arrayLoessCurveXY(1,:), arrayLoessCurveXY(2,:), arrayOutXCentresRescaled, 'pchip');
        

end

