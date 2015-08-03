function [ arrAbundVect ] = convertSmpAnaToAbundVector( structSampleInfo, numSpatialBins, arrayNormDistRatio, numLoessWindowSize, flagPerformZNorm )
%convertSmpAnaToAbundVector
%   Detailed explanation goes here

    [ arrayOutXCentres, ~, ~ ] = calculateSpatialDivisions( numSpatialBins, arrayNormDistRatio );
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

    arraySampleMean = mean(double(arrayTempSampleValues(:)));
    arraySampleStDev = std(double(arrayTempSampleValues(:)));

    if flagPerformZNorm,
        arrayDataToSmooth = (double(arrayTempSampleValues) - arraySampleMean)/arraySampleStDev;
    else
        arrayDataToSmooth = arrayTempSampleValues;
    end

    arrayRescaledSamplePositions = rescaleNormDist( arrayTempSamplePositions, numSpatialBins, arrayNormDistRatio );

    arrayLoessCurveXY = calculateLoessCurve( arrayRescaledSamplePositions, arrayDataToSmooth, numSpatialBins, arrayNormDistRatio, numLoessWindowSize );
    
    arrAbundVect = interp1(arrayLoessCurveXY(1,:), arrayLoessCurveXY(2,:), arrayOutXCentresRescaled, 'pchip');
        

end

