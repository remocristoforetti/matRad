function [outputProbabilisticQuantities,probQuantitiesAccumulationTime] = matRad_accumulateProbabilisticQuantities(saveDir,ct,cst,scenariosMeta, multiScen, mode4D, bioQuantities)
    
    matRad_cfg = MatRad_Config.instance();
    originaLogLevel = matRad_cfg.logLevel;
    testVariance = false;

    % Input handle
    if ~exist('mode4D', 'var') || isempty(mode4D)
        mode4D = 'all';
    elseif ~any(strcmpi(mode4D,{'phase','all'}))
        matRad_cfg.dispError('4D calculation mode %s unknown, allowed is only ''phase'' or ''all''',mode4D);
    end

    if ~exist('scenariosMeta', 'var')
        scenariosMeta = matRad_getMetaFromScenariosPool(saveDir);
    end

    
    if ~exist('multiScen', 'var')
        matRad_cfg.dispWarning('Specific multScen not provided, loading default one');
        multiScen = matRad_getMultiScenFromScenarios(saveDir,'rndScen');
    end

    if ~exist('bioQuantities', 'var') || isempty(bioQuantities)
        bioQuantities = false;
    end

    % Load dij template
    try

        load(fullfile(saveDir, 'dijTemplate'));
    catch
        matRad_cfg.dispError('Unable to load dijTemplate file');
    end

    % Check consistency
    nScenarios = numel(scenariosMeta);

    if nScenarios ~= multiScen.totNumScen
        matRad_cfg.dispError('Inconsistent number of scenarios detected');
    end

    
    % mode selection
    switch mode4D
        case 'phase'
            for ctIx = 1:multiScen.numOfCtScen

                scenIx = multiScen.linearMask(:,1) == ctIx;
                
                ctAccumIx{ctIx}(:,1:3) = multiScen.linearMask(scenIx,:);
                
                % This is an index corresponding to the scenario position
                % in the scenarioMeta
                ctAccumIx{ctIx}(:,4) = find(scenIx);

                % get weights associated to thes scenarios;
                scenWeights{ctIx} = multiScen.scenWeight(scenIx);
                
                %normalize weights for current ctScenario
                scenWeights{ctIx} = scenWeights{ctIx}./sum(scenWeights{ctIx});
   
            end
                
        case 'all'
            % This assumes that multScen has a
            ctAccumIx{1} = multiScen.linearMask;
            ctAccumIx{1}(:,4) = [1:nScenarios];
            scenWeights{1} = multiScen.scenWeight./sum(multiScen.scenWeight);
    
    end

    structsToInclude = find(~cellfun(@isempty, cst(:,6)))';
    cst = matRad_setOverlapPriorities(cst,dijTemplate.doseGrid.dimensions);

    % Memory Allocation
    tic;
    switch mode4D
        case 'phase'
            
            expDist = cell(multiScen.numOfCtScen,1);
            expDist(:) = {spalloc(dijTemplate.doseGrid.numOfVoxels,dijTemplate.totalNumOfBixels,1)};
            omega = cell(size(cst,1),multiScen.numOfCtScen);
            omega(:) = {spalloc(dijTemplate.totalNumOfBixels, dijTemplate.totalNumOfBixels,1)};

            if bioQuantities
                mAlphaDoseExp    = cell(multiScen.numOfCtScen,1);
                mAlphaDoseExp(:) = {spalloc(dijTemplate.doseGrid.numOfVoxels,dijTemplate.totalNumOfBixels,1)};
    
                mSqrtBetaDoseExp    = cell(multiScen.numOfCtScen,1);
                mSqrtBetaDoseExp(:) = {spalloc(dijTemplate.doseGrid.numOfVoxels,dijTemplate.totalNumOfBixels,1)};
    
                mAlphaDoseOmega    = cell(size(cst,1),multiScen.numOfCtScen);
                mAlphaDoseOmega(:) = {spalloc(dijTemplate.totalNumOfBixels, dijTemplate.totalNumOfBixels,1)};
    
                mSqrtBetaDoseOmega    = cell(size(cst,1),multiScen.numOfCtScen);
                mSqrtBetaDoseOmega(:) = {spalloc(dijTemplate.totalNumOfBixels, dijTemplate.totalNumOfBixels,1)};

            end

  
        case 'all'
            expDist = {spalloc(dijTemplate.doseGrid.numOfVoxels,dijTemplate.totalNumOfBixels,1)};
            omega   = cell(size(cst,1),1);
            omega(:) = {spalloc(dijTemplate.totalNumOfBixels, dijTemplate.totalNumOfBixels,1)};
            
            physicalDoseJExp = cell(size(cst,1),1);
            physicalDoseJExp(:) = {zeros(dijTemplate.totalNumOfBixels,1)};
            
            if bioQuantities
                mAlphaDoseExp    = {spalloc(dijTemplate.doseGrid.numOfVoxels,dijTemplate.totalNumOfBixels,1)};
                mSqrtBetaDoseExp = {spalloc(dijTemplate.doseGrid.numOfVoxels,dijTemplate.totalNumOfBixels,1)};
    
                mAlphaDoseOmegaExp   = cell(size(cst,1),1);
                mAlphaDoseOmegaExp(:) = {spalloc(dijTemplate.totalNumOfBixels, dijTemplate.totalNumOfBixels,1)};
                 
                mSqrtBetaDoseOmegaExp   = cell(size(cst,1),1);
                mSqrtBetaDoseOmegaExp(:) = {spalloc(dijTemplate.totalNumOfBixels, dijTemplate.totalNumOfBixels,1)};
    
                % mAlphaSqrtBetaDoseOmega = cell(size(cst,1),1);
                % mAlphaSqrtBetaDoseOmega(:) = {spalloc(dijTemplate.totalNumOfBixels, dijTemplate.totalNumOfBixels,1)};
                % 
                % mTwiceAlphaSqrtBetaDoseOmega = cell(size(cst,1),1);
                % mTwiceAlphaSqrtBetaDoseOmega(:) = {spalloc(dijTemplate.totalNumOfBixels, dijTemplate.totalNumOfBixels,1)};
    
                alphaJExp = cell(size(cst,1),1);
                alphaJExp(:) = {zeros(dijTemplate.totalNumOfBixels,1)};
    
                sqrtBetaJExp = cell(size(cst,1),1);
                sqrtBetaJExp(:) = {zeros(dijTemplate.totalNumOfBixels,1)};

                % mAlphaDoseOmegaExp   = cell(size(cst,1),1);
                % mAlphaDoseOmegaExp(:) = {spalloc(dijTemplate.totalNumOfBixels, dijTemplate.totalNumOfBixels,1)};

                % mAlphaDoseOmegaCross   = cell(size(cst,1),1);
                % mAlphaDoseOmegaCross(:) = {spalloc(dijTemplate.totalNumOfBixels, dijTemplate.totalNumOfBixels,1)};

                % mSqrtBetaDoseOmegaCross   = cell(size(cst,1),1);
                % mSqrtBetaDoseOmegaCross(:) = {spalloc(dijTemplate.totalNumOfBixels, dijTemplate.totalNumOfBixels,1)};
            end
    end


    gpu = gpuDevice();
    
    matRad_cfg.dispInfo('Accumulating probabilistic quantities E[D] & Omega[D] ...\n');


    for phaseIdx=1:numel(ctAccumIx)
        lineLengthPhase = fprintf('\t\t4D-Phase %d/%d...\n',phaseIdx,numel(ctAccumIx));

        currPhaseOmega = omega(:,phaseIdx);
        currPhasePhysicalDoseJ = physicalDoseJExp(:,phaseIdx);

        if bioQuantities

            currPhaseOmegaAlpha = mAlphaDoseOmegaExp(:,phaseIdx);
            currPhaseOmegaBeta  = mSqrtBetaDoseOmegaExp(:,phaseIdx);
            
            % currPhaseOmegaAlphaCross = mAlphaDoseOmegaCross(:,phaseIdx);
            % currPhaseOmegaSqrtBetaDose = mSqrtBetaDoseOmegaCross(:,phaseIdx);
            
            currPhaseAlphaJ     = alphaJExp(:,phaseIdx);
            currPhaseSqrtBetaJ  = sqrtBetaJExp(:,phaseIdx);
            


        end

        % Get scenarios meta in this psecific phase
        scenariosMetaInPhase = scenariosMeta(ctAccumIx{phaseIdx}(:,4));
        nScenariosInPhase = numel(scenariosMetaInPhase);

        % Accumulate quantities
        for scenIdx=1:nScenariosInPhase
            lineLengthScen = fprintf('\t\t\tAccumulating scenario: %d/%d\n',scenIdx,nScenariosInPhase);
            currMeta = scenariosMetaInPhase(scenIdx);
            outScenarioQuantities = matRad_loadScenariosFromMeta(saveDir,currMeta,dijTemplate, bioQuantities);

            scenarioDistribution = outScenarioQuantities.physicalDose;
            scenarioPhysicalDoseJ = cell(size(cst,1),1);
            scenarioPhysicalDoseJ(structsToInclude)       = arrayfun(@(structIdx) sum(scenarioDistribution{1}(cst{structIdx,4}{currMeta.ctScenIdx},:),1), structsToInclude, 'UniformOutput',false);

            lineLengthScen = lineLengthScen + fprintf('\t\t\t\tAccumulating exp...');
            expDist{phaseIdx} = expDist{phaseIdx} + scenarioDistribution{1}.*scenWeights{phaseIdx}(scenIdx);

            if bioQuantities
  
                ax = zeros(size(scenarioDistribution{1},1),1);
                bx = zeros(size(scenarioDistribution{1},1),1);
    
                for v = 1:size(cst,1)
                    ax(cst{v,4}{currMeta.ctScenIdx}) = cst{v,5}.alphaX;
                    bx(cst{v,4}{currMeta.ctScenIdx}) = cst{v,5}.betaX;
                end

                if isfield(outScenarioQuantities, 'mAlphaDose')  && nnz(outScenarioQuantities.mAlphaDose{1})>0
                    scenarioAlphaDose    = outScenarioQuantities.mAlphaDose;
                else
                    scenarioAlphaDose    = {scenarioDistribution{1}.*ax};
                end
                
                if isfield(outScenarioQuantities, 'mSqrtBetaDose') && nnz(outScenarioQuantities.mSqrtBetaDose{1})>0
                    scenarioSqrtBetaDose = outScenarioQuantities.mSqrtBetaDose;
                else
                    scenarioSqrtBetaDose    = {scenarioDistribution{1}.*sqrt(bx)};
                end

                scenarioAlphaJ = cell(size(cst,1),1);
                scenarioSqrtBetaJ = cell(size(cst,1),1);
                
                scenarioAlphaJ(structsToInclude)       = arrayfun(@(structIdx) sum(scenarioAlphaDose{1}(cst{structIdx,4}{currMeta.ctScenIdx},:),1), structsToInclude, 'UniformOutput',false);
                scenarioSqrtBetaJ(structsToInclude)    = arrayfun(@(structIdx) sum(scenarioSqrtBetaDose{1}(cst{structIdx, 4}{currMeta.ctScenIdx},:),1), structsToInclude, 'UniformOutput',false);
                
                mAlphaDoseExp{phaseIdx}    = mAlphaDoseExp{phaseIdx} + scenarioAlphaDose{1}.*scenWeights{phaseIdx}(scenIdx);
                mSqrtBetaDoseExp{phaseIdx} = mSqrtBetaDoseExp{phaseIdx} + scenarioSqrtBetaDose{1}.*scenWeights{phaseIdx}(scenIdx);

            end

            lineLengthScen = lineLengthScen + fprintf('done.\n');
            lineLengthScen = lineLengthScen + fprintf('\t\t\t\tAccumulating omega for struct: ');
            
            structCounter = 0;
            
            for structIdx=structsToInclude
                structCounter = structCounter+1;
                lineLengthStruct = fprintf('%d/%d',structCounter,numel(structsToInclude));

                % This is different from previous branch. I consider here
                % only voxels within the structure in this ct scenario.
                % Beore I was considering dose to structure in all ct
                % scenarios
                %currStructVoxels = cst{structIdx,4}{currMeta.ctScenIdx};

                % Old branch
                % Get all voxels in current structure on all ct scenarios
                scenariosMetaInPhase = scenariosMeta(ctAccumIx{phaseIdx}(:,4));
                nScenariosInPhase = numel(scenariosMetaInPhase);
                currStructVoxels = [];
                for scenIdx=1:nScenariosInPhase
                   currMeta = scenariosMetaInPhase(scenIdx);
                   currStructVoxels = [currStructVoxels; cst{structIdx,4}{currMeta.ctScenIdx}];
                end
                currStructVoxels = unique(currStructVoxels);

                currStructOmega = currPhaseOmega{structIdx};
                currDist = gpuArray(scenarioDistribution{1}(currStructVoxels,:));

                currStructOmega = currStructOmega + gather(currDist'*currDist).*scenWeights{phaseIdx}(scenIdx);
                wait(gpu);

                currPhaseOmega{structIdx} = currStructOmega;
                
                currStructPhysicalDoseJ = currPhasePhysicalDoseJ{structIdx};
                currStructPhysicalDoseJ = currStructPhysicalDoseJ + scenarioPhysicalDoseJ{structIdx}'.*scenWeights{phaseIdx}(scenIdx);                    
                currPhasePhysicalDoseJ{structIdx} = currStructPhysicalDoseJ;

                if bioQuantities
                    % omega Alpha
                    currStructOmegaAlpha = currPhaseOmegaAlpha{structIdx};
                    currDistAlpha = gpuArray(scenarioAlphaDose{1}(currStructVoxels,:));
    
                    currStructOmegaAlpha = currStructOmegaAlpha + gather(currDistAlpha'*currDistAlpha).*scenWeights{phaseIdx}(scenIdx);
                    wait(gpu);
                    
                    currPhaseOmegaAlpha{structIdx} = currStructOmegaAlpha;
    
                    % Omega Beta
                    currStructOmegaBeta = currPhaseOmegaBeta{structIdx};
                    currDistBeta = gpuArray(scenarioSqrtBetaDose{1}(currStructVoxels,:));
                    currStructOmegaBeta = currStructOmegaBeta + gather(currDistBeta'*currDistBeta).*scenWeights{phaseIdx}(scenIdx);
                    wait(gpu);
    
                    currPhaseOmegaBeta{structIdx} = currStructOmegaBeta;

                    % AlphaJ
                    currStructAlphaJ = currPhaseAlphaJ{structIdx};
                    currStructAlphaJ = currStructAlphaJ + scenarioAlphaJ{structIdx}'.*scenWeights{phaseIdx}(scenIdx); % This is the scenario average
                    
                    currPhaseAlphaJ{structIdx} = currStructAlphaJ;

                    % % AlphaCross
                    % currStructOmegaAlphaCross = currPhaseOmegaAlphaCross{structIdx};
                    % currAlphaOmegaCross =  scenarioAlphaJ{structIdx}'*scenarioAlphaJ{structIdx};
                    % 
                    % currStructOmegaAlphaCross = currStructOmegaAlphaCross + currAlphaOmegaCross.*scenWeights{phaseIdx}(scenIdx);
                    % wait(gpu);
                    % currPhaseOmegaAlphaCross{structIdx} = currStructOmegaAlphaCross;

                    % % Sqrt Beta Croos
                    % currStructOmegaSqrtBetaCross = currPhaseOmegaSqrtBetaDose{structIdx};
                    % currSqrtBetaOmegaCross =  scenarioSqrtBetaJ{structIdx}'*scenarioSqrtBetaJ{structIdx};
                    % 
                    % currStructOmegaSqrtBetaCross = currStructOmegaSqrtBetaCross + currSqrtBetaOmegaCross.*scenWeights{phaseIdx}(scenIdx);
                    % wait(gpu);
                    % currPhaseOmegaSqrtBetaDose{structIdx} = currStructOmegaSqrtBetaCross;

                    % BetaJ
                    currStructSqrtBetaJ = currPhaseSqrtBetaJ{structIdx};
                    currStructSqrtBetaJ = currStructSqrtBetaJ + scenarioSqrtBetaJ{structIdx}'*scenWeights{phaseIdx}(scenIdx);
    
                    currPhaseSqrtBetaJ{structIdx} = currStructSqrtBetaJ;
                end

                clear currDist;
                clear currDistAlpha;
                clear currDistBeta;

                if structIdx~=structsToInclude(end)
                    fprintf(repmat('\b',1,lineLengthStruct));
                    lineLengthStruct = 0;
                end
                
            end
            lineLengthScen = lineLengthScen + fprintf(' done.\n');

            if scenIdx~=nScenariosInPhase
                fprintf(repmat('\b',1,lineLengthScen+4));
                fprintf('\n');
                lineLengthScen = 0;
            end
        end

        outputProbabilisticQuantities.physicalDoseOmegaExp(:,phaseIdx)            = currPhaseOmega;
        outputProbabilisticQuantities.physicalDoseJExp(:,phaseIdx)                = currPhasePhysicalDoseJ;
       
        if bioQuantities

            outputProbabilisticQuantities.mAlphaDoseOmegaExp(:,phaseIdx)         = currPhaseOmegaAlpha;
            outputProbabilisticQuantities.mSqrtBetaDoseOmegaExp(:,phaseIdx)      = currPhaseOmegaBeta;
            outputProbabilisticQuantities.alphaDoseJExp(:, phaseIdx)             = currPhaseAlphaJ;
            outputProbabilisticQuantities.sqrtBetaDoseJExp(:,phaseIdx)           = currPhaseSqrtBetaJ;


            % for structIdx=structsToInclude
            % 
            %     % Get voxel indexes
            %     scenariosMetaInPhase = scenariosMeta(ctAccumIx{phaseIdx}(:,4));
            %     nScenariosInPhase = numel(scenariosMetaInPhase);
            %     currStructVoxels = [];
            %     for scenIdx=1:nScenariosInPhase
            %        currMeta = scenariosMetaInPhase(scenIdx);
            %        currStructVoxels = [currStructVoxels; cst{structIdx,4}{currMeta.ctScenIdx}];
            %     end
            %     currStructVoxels = unique(currStructVoxels);
            % 
            %     currDistAlphaExp = gpuArray(mAlphaDoseExp{1}(currStructVoxels,:));
            %     mAlphaDoseOmegaExp{structIdx} = gather(currDistAlphaExp'*currDistAlphaExp); % Sum over voxels of expAlphaDose
            %     wait(gpu);
            % 
            % end
            %     outputProbabilisticQuantities.mAlphaDoseOmegaExp(:,phaseIdx) = mAlphaDoseOmegaExp;
            %     outputProbabilisticQuantities.mAlphaDoseOmegaCross(:,phaseIdx) = currPhaseOmegaAlphaCross;
            %     outputProbabilisticQuantities.mSqrtBetaDoseOmegaCross(:,phaseIdx) = currPhaseOmegaSqrtBetaDose;        
        end
        
        if phaseIdx~=numel(ctAccumIx)
            fprintf(repmat('\b',1,lineLengthPhase));
            lineLengthPhase = 0;
        end

    end

    outputProbabilisticQuantities.physicalDoseExp  = expDist;

    if bioQuantities
        outputProbabilisticQuantities.mAlphaDoseExp    = mAlphaDoseExp;
        outputProbabilisticQuantities.mSqrtBetaDoseExp = mSqrtBetaDoseExp;
        
    end

    matRad_cfg.dispInfo('done.\n');
    probQuantitiesAccumulationTime = toc;

end