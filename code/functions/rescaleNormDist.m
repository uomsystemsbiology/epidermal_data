function [ arrayOutNormDistRescaled ] = rescaleNormDist( arrayInNormDist, numInSpatialBins, arrayInNormDistRatio )
%rescaleNormDist :: Take the input array (normalised distance values) and
%rescale them for plotting according to to the number of spatial bins

    numSamples = size(arrayInNormDist,1);
    numSamplesPerRatio = numInSpatialBins/sum(arrayInNormDistRatio);
    
    arrayOutNormDistRescaled = double(zeros(numSamples,1));
    
    for iSample = 1:numSamples,
        %identify the layer from the whole-number component
        numSampleNormDist = arrayInNormDist(iSample);
        numSampleLayer = floor(numSampleNormDist);
        
        %the fractional distance is then scaled according to the layer
        %distance
        numSampleFractionalDistance = numSampleNormDist-numSampleLayer;
        if numSampleLayer == 0,
            numSampleNormDistFractionRescaled = numSampleFractionalDistance*numSamplesPerRatio*arrayInNormDistRatio(1);
            numSampleNormDistLayerRescaled = 0;
        elseif numSampleLayer == 1,
            numSampleNormDistFractionRescaled = numSampleFractionalDistance*numSamplesPerRatio*arrayInNormDistRatio(2);
            numSampleNormDistLayerRescaled = numSamplesPerRatio*arrayInNormDistRatio(1);
        elseif numSampleLayer == 2,
            numSampleNormDistFractionRescaled = numSampleFractionalDistance*numSamplesPerRatio*arrayInNormDistRatio(3);
            numSampleNormDistLayerRescaled = numSamplesPerRatio*sum(arrayInNormDistRatio(1:2));
        elseif numSampleLayer == 3,
            numSampleNormDistFractionRescaled = 0;
            numSampleNormDistLayerRescaled = numSamplesPerRatio*sum(arrayInNormDistRatio(1:3));
        end
        
        numSampleNormDistRescaled = numSampleNormDistLayerRescaled + numSampleNormDistFractionRescaled;
        
        arrayOutNormDistRescaled(iSample) = numSampleNormDistRescaled;
    end

end

