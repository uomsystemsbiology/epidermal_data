%% ==============================================================
% Define Input Variables
% ===============================================================

% ------- Define MATLAB paths -------
stringSharedScriptsFolder = 'C:\wc\code\shared_functions';
addpath(genpath(stringSharedScriptsFolder));

stringCurrentDirectory = cd;

%check that the network inference framework has been initialised
global SessionState;
if ~exist('SessionState.frameworkRoot', 'var'),
    cd(getenv('$framework_root'));
    run scripts\setup;
    cd(stringCurrentDirectory);
end
SessionState.isRunningStandalone = true;
% -----------------------------------


% ------ Define the input data domain ------
numTotalNodes = 51;
arrayOutputPatients = 1:3;
numPatients = length(arrayOutputPatients);

arrayOutputNodes = 1:15; %[ 1:17 ];
numOutputNodes = length(arrayOutputNodes);

numTissueLayers = 4;

%read in an empty network (combined format, has all nodes)
pathEmptyNetworkFile = [ stringCurrentDirectory '\emptyNetwork.txt'];
arrayEmptyNetwork = readNetwork(pathEmptyNetworkFile,'-combined');

% -----------------------------------

numTOSTAlpha = 0.05;
numTOSTBeta = numTOSTAlpha;
numTOSTDelta = 0;

% ------ Define the boundary file names ------
arrayBoundFileStrings = { [ 'im_basal_memb_z' ];
                          [ 'im_bound1_z' ];
                          [ 'im_bound2_z' ];
                          [ 'im_outer_bound_z' ] };
numLayerBounds = length(arrayBoundFileStrings);
arrayLayerBoundCombinations = { [1 2];
                                [2 3];
                                [3 4];
                                [1 4] };
arrayLayerStrings = { {['Basal'];['Layer']};
                      {['Spinous and'];['Granular Layers']};
                      {['Transitional'];['Layer']};
                      {['Total'];['Epidermis']} };
                  
arraySubFigStrings = { ['(a)'];
                       ['(b)'];
                       ['(c)'];
                       ['(d)'] };
% -----------------------------------


arrayPatientColors = { [ 1 0 0 ];
                       [ 0 1 0 ];
                       [ 0 0 1 ] };


% ------ Define the image pixel resolutions (um per pixel) ------
arrayPixRes = [ 0.148160, 0.157981, 0.144300;       %ITGB1
                0.131244, 0.161216, 0.151396;       %EGFR_ph
                0.142824, 0.154177, 0.170356;       %ITGB4
                0.152020, 0.144811, 0.160251;       %SFN
                0.135444, 0.132890, 0.128008;       %CALM
                0.168028, 0.173308, 0.154915;       %RAF
                0.105302, 0.137147, 0.162636;       %RAF_ph
                0.139532, 0.134139, 0.119153;       %MEK
                0.142200, 0.172343, 0.120401;       %MEK_ph
                0.147479, 0.159570, 0.149920;       %ERK
                0.145719, 0.232515, 0.151396;       %ERK_ph
                0.140951, 0.141575, 0.129257;       %JUNB
                0.138680, 0.154121, 0.129257;       %JUNC
                0.149750, 0.156051, 0.110638;       %FOSC
                0.144868, 0.163601, 0.118187;       %FRA2
                0.115122, 0.132266, 0.149012;       %K10
                0.125964, 0.137829, 0.126135 ];     %K14

% -----------------------------------

%% ==============================================================
% Perform Pre-Processing
% ===============================================================

%check that the folder names (stored in 'arrayOutputFolders') are initialised
if ~exist('arrayOutputFolders', 'var'),
    [ arrayOutputFolders ] = initialise_output_folders(  );
end



%% ==============================================================
% Load the layer-boundary data and minimum distances
% ===============================================================
arrayLayerPixThickness = cell(numPatients, numOutputNodes, numTissueLayers);
arrayZSlices = cell(numPatients,numOutputNodes);
for iPatient = 1:numPatients,
    numPatient = arrayOutputPatients(iPatient);
    disp(['Patient ' num2str(numPatient)]);
    
    for iNode = 1:numOutputNodes,
        numNode = arrayOutputNodes(iNode);
        disp(['Node ' num2str(numNode)]);
        
        stringDataFolder = arrayOutputFolders{numNode, numPatient};
        arrayBackslashLoc = strfind(stringDataFolder, '\');
        stringImageFolder = stringDataFolder(1:(arrayBackslashLoc(end-1)));
        
        %load the new data and extract the z-positions
        clear SampleAnalysis;
        load(  [ stringDataFolder 'sample_ana.mat' ]  );
        arrayZSlices{iPatient,iNode} = find_z_slices_sampled( SampleAnalysis );
        numZSlices = length(arrayZSlices{iPatient,iNode});
        for iLayer = 1:numTissueLayers,
            arrayLayerPixThickness{iPatient, iNode, iLayer} = cell(numZSlices);
        end
        
        %move through each sampled z-position
        for iZPos = 1:numZSlices,
        
            numZPos = arrayZSlices{iPatient,iNode}(iZPos);
            
            %move through each tissue layer to measure the minimum distance
            for iLayer = 1:numTissueLayers,
                
                
                %load the masked boundary images
                imageLowerBound = imread( [stringImageFolder arrayBoundFileStrings{arrayLayerBoundCombinations{iLayer}(1)} num2str(numZPos) '.tiff']  );
                imageUpperBound = imread( [stringImageFolder arrayBoundFileStrings{arrayLayerBoundCombinations{iLayer}(2)} num2str(numZPos) '.tiff']  );
                
                %extract the pixel co-ordinates
                [arrayLowerY, arrayLowerX] = find(imageLowerBound);
                [arrayUpperY, arrayUpperX] = find(imageUpperBound);
                numLowerPix = length(arrayLowerY);
                numUpperPix = length(arrayUpperY);
                
                %create the output vector
                arrayLayerPixThickness{iPatient, iNode, iLayer}{iZPos} = zeros(numLowerPix,1,'double');
                
                %combine the pixel lists and use the dist function
                arrayCombinedPix = [ arrayLowerX', arrayUpperX'; arrayLowerY', arrayUpperY' ];
                arrayFullDistMatrix = dist(arrayCombinedPix);
                
                %extract the section of the distance matrix which is of interest
                arrayLowerToUpper = arrayFullDistMatrix(1:numLowerPix,numLowerPix+1:end);
                
                %and determine the minimum values
                for iPix = 1:numLowerPix,
                    arrayLayerPixThickness{iPatient, iNode, iLayer}{iZPos}(iPix) = min(arrayLowerToUpper(iPix,:));
                end
                
            end
            
        end
        
    end
    
end


figCompPatients = figure;
hold on;
arrayCombAllObs = cell(numPatients,numTissueLayers);
for iLayer = 1:numTissueLayers,
    
    subplot(2,2,iLayer);
    hold on;
    
    numCombObs = zeros(numPatients,1,'double');
    numCalcMean = zeros(numPatients,1,'double');
    numCalcStDev = zeros(numPatients,1,'double');
    
    for iPatient = 1:numPatients,
        
        %create the combined/aggregate array
        numObs = 0;
        for iNode = 1:numOutputNodes,
            for iZPos = 1:length(arrayZSlices{iPatient,iNode}),
               numObs = numObs + length( arrayLayerPixThickness{iPatient, iNode, iLayer}{iZPos} );
            end
        end
        arrayCombAllObs{iPatient,iLayer} = zeros(numObs,1,'double');
        
        %populate the combined/aggregate array
        numStart = 1;
        for iNode = 1:numOutputNodes,
            for iZPos = 1:length(arrayZSlices{iPatient,iNode}),
               numObs = length( arrayLayerPixThickness{iPatient, iNode, iLayer}{iZPos} );
               numEnd = numStart+numObs-1;
               arrayCombAllObs{iPatient,iLayer}(numStart:numEnd) = arrayLayerPixThickness{iPatient, iNode, iLayer}{iZPos};
               numStart = numEnd+1;
            end
        end
        
        %convert from pixel distance to absolute distance
        arrayCombAllObs{iPatient,iLayer} = arrayCombAllObs{iPatient,iLayer}*arrayPixRes(iNode,numPatient);
        [arrayProbDens, arrayDist] = ksdensity(arrayCombAllObs{iPatient,iLayer});
        
        numCombObs(iPatient) = length(arrayCombAllObs{iPatient,iLayer});
        numCalcMean(iPatient) = mean(arrayCombAllObs{iPatient,iLayer});
        numCalcStDev(iPatient) = std(arrayCombAllObs{iPatient,iLayer});
        
        plot(arrayDist, arrayProbDens, '-', 'Color', arrayPatientColors{iPatient}, 'LineWidth', 2);
        
    end
    arrayXLim = get(gca, 'XLim');
    set(gca, 'XLim', [0 arrayXLim(2)]);
    arrayYLim = get(gca, 'YLim');
    set(gca, 'YLim', [0 arrayYLim(2)]);
    set(gca, 'FontSize', 14, 'FontName', 'Times New Roman');
    text(0.65*arrayXLim(2),0.85*arrayYLim(2),arrayLayerStrings{iLayer}, 'FontSize', 14, 'FontName', 'Times New Roman', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    text(-0.3*arrayXLim(2), 1.05*arrayYLim(2), arraySubFigStrings{iLayer}, 'FontSize', 14, 'FontName', 'Times New Roman');
    xlabel('Abs. Dist. ({\mu}m)', 'FontSize', 14, 'FontName', 'Times New Roman');
    ylabel('pdf', 'FontSize', 14, 'FontName', 'Times New Roman');
    
   
end
hold off;
print(figCompPatients, '-r300', '-dpng', [ stringCurrentDirectory '\PatALL_-_CompThickness.png' ]);
close(figCompPatients);



arrayTempCols = jet(numOutputNodes);

arrayCombNodeObs = cell(numPatients,numOutputNodes,numTissueLayers);
for iPatient = 1:numPatients,
    
    figWithinPat = figure;
	hold on;
    
    for iLayer = 1:numTissueLayers,
    
        subplot(2,2,iLayer);
        hold on;
    
        for iNode = 1:numOutputNodes,
        
            %combine the data over all z-layers
            numObs = 0;
            for iZPos = 1:length(arrayZSlices{iPatient,iNode}),
               numObs = numObs + length( arrayLayerPixThickness{iPatient, iNode, iLayer}{iZPos} );
            end
            arrayCombNodeObs{iPatient,iNode,iLayer} = zeros(numObs,1,'double');
            
        
            %populate the combined/aggregate array
            numStart = 1;
            for iZPos = 1:length(arrayZSlices{iPatient,iNode}),
               numObs = length( arrayLayerPixThickness{iPatient, iNode, iLayer}{iZPos} );
               numEnd = numStart+numObs-1;
               arrayCombNodeObs{iPatient,iNode,iLayer}(numStart:numEnd) = arrayLayerPixThickness{iPatient, iNode, iLayer}{iZPos};
               numStart = numEnd+1;
            end
        
            %convert from pixel distance to absolute distance
            arrayCombNodeObs{iPatient,iNode,iLayer} = arrayCombNodeObs{iPatient,iNode,iLayer}*arrayPixRes(iNode,numPatient);
            [arrayProbDens, arrayDist] = ksdensity(arrayCombNodeObs{iPatient,iNode,iLayer});

            plot(arrayDist, arrayProbDens, '-', 'Color', arrayTempCols(iNode,:));
        
        end
        
        title(arrayLayerStrings{iLayer}, 'FontName', 'Times New Roman');
        xlabel('Abs. Dist. (um)', 'FontName', 'Times New Roman');
        ylabel('pdf', 'FontName', 'Times New Roman');
        arrayXLim = get(gca, 'XLim');
        set(gca, 'XLim', [0 arrayXLim(2)]);
        
        %compare different nodes within each tissue layer
        
        
    
    end

    
    print(figWithinPat, '-r300', '-dpng', [stringCurrentDirectory '\Pat' num2str(iPatient) '_CompThickness.png']);
    close(figWithinPat);
    
end

