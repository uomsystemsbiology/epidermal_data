function flagResult = extractSampleInformation(structSettings)

    %extract_sample_information.m
    %
    %   This script moves through the specified target proteins/sub-cellular
    %localisations and populates the SampleAnalysis structured arrays with
    %information on sample location, associated signal intensity, and
    %spatial-normalisation (i.e. the calculated layer-normalised distance)
    %output:
    %   - SampleAnalysis: a structured array which is used for subsequent
    %           data analysis, containing the fields: 
    %           - CoOrds: a vector with the x, y and z co-ordinates of the
    %                   sampled pixel
    %           - Layer: an integer with the tissue layer
    %           - BasDist: the distance (in pixels) from the basement membrane
    %           - NormDist: layer-normalized distance (linear interpolation
    %                   within each tissue layer)
    %           - SigInt: the signal intensity (uint8) associated with the
    %                   pixel




    %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
    %% Produce the Output Arrays
    %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
    %move through each target protein/phospho-protein
    for iTarget = 1:numTargets,

        %move through each patient
        for iPatient = 1:numPatients,

            %move through each sub-cellular localisation
            for iLoc = 1:numLocalisations
                %load the sample_loc.mat array
                stringSegImgPath = [ strProcImgDataPath arrayTargetPaths{iTarget} '\Pat_' num2str(iPatient) '\' ];
                stringSmpLocPath = [ strProcImgDataPath arrayTargetPaths{iTarget} '\Pat_' num2str(iPatient) '\' arrayLocStrings{iLoc} '\' ];
                load([stringSmpLocPath 'sample_loc.mat']);
                structSampleLocs = SampleOutput;
                clear SampleOutput;

                %determine the path of the image stacks
                stringImgDataPath = [ strRawImgDataPath arrayTargetPaths{iTarget} '\Pat_' num2str(iPatient) '\image_data\' ];

                %produce a 'SampleAnalysis' structured array
                structSignalIntensityData = produceSmpAna(structSampleLocs, arraySamplingKernel, stringImgDataPath, stringSegImgPath);

                %save this structured array
                SampleAnalysis = structSignalIntensityData;
                save([stringSmpLocPath 'sample_ana.mat'], 'SampleAnalysis');
                clear SampleAnalysis;

            end

        end


    end
    
    flagResult = true;
    
end