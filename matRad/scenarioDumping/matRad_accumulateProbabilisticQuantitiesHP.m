function [outputProbabilisticQuantities,probQuantitiesAccumulationTime] = matRad_accumulateProbabilisticQuantitiesHP(saveDir,ct,cst,scenariosMeta, multiScen, mode4D, bioQuantities)
    
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
    
                % This is alpha(ij) * alpha(ik)
                mAlphaDoseOmegaExp   = cell(size(cst,1),1);
                mAlphaDoseOmegaExp(:) = {spalloc(dijTemplate.totalNumOfBixels, dijTemplate.totalNumOfBixels,1)};
                
                % This is sqrtBeta(ij) * sqrtBeta(ik)
                mSqrtBetaDoseOmegaExp   = cell(size(cst,1),1);
                mSqrtBetaDoseOmegaExp(:) = {spalloc(dijTemplate.totalNumOfBixels, dijTemplate.totalNumOfBixels,1)};
                
                
                % This is alpha(ij) * alpha(ik) + 2*alpha(ij)*beta(ik) + beta(ij)*beta(ik)
                % -> Beta here is square of sqrtBeta. Only trace of  the sqrtBeta x sqrtBeta ,atrix, negleting covariance terms in sqrtBeta 
                mEffectOmegaHP = cell(size(cst,1),1);
                mEffectOmegaHP(:) = {spalloc(dijTemplate.totalNumOfBixels, dijTemplate.totalNumOfBixels,1)};
            end
    end


    gpu = gpuDevice();
    
    matRad_cfg.dispInfo('Accumulating probabilistic quantities E[D] & Omega[D] ...\n');


    for phaseIdx=1:numel(ctAccumIx)
        lineLengthPhase = fprintf('\t\t4D-Phase %d/%d...\n',phaseIdx,numel(ctAccumIx));

        currPhaseOmega = omega(:,phaseIdx);


        if bioQuantities

            currPhaseOmegaAlpha = mAlphaDoseOmegaExp(:,phaseIdx);
            currPhaseOmegaBeta  = mSqrtBetaDoseOmegaExp(:,phaseIdx);
            
            currPhaseOmegaEffect = mEffectOmegaHP(:,phaseIdx);

        end

        % Get scenarios meta in this psecific phase
        scenariosMetaInPhase = scenariosMeta(ctAccumIx{phaseIdx}(:,4));
        nScenariosInPhase    = numel(scenariosMetaInPhase);

        % Accumulate quantities
        for scenIdx=1:nScenariosInPhase
            lineLengthScen = fprintf('\t\t\tAccumulating scenario: %d/%d\n',scenIdx,nScenariosInPhase);
            currMeta = scenariosMetaInPhase(scenIdx);
            outScenarioQuantities = matRad_loadScenariosFromMeta(saveDir,currMeta,dijTemplate, bioQuantities);

            scenarioDistribution = outScenarioQuantities.physicalDose;

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

                scenarioBetaDose = {scenarioSqrtBetaDose{1}.^2};
                
                mAlphaDoseExp{phaseIdx}    = mAlphaDoseExp{phaseIdx} + scenarioAlphaDose{1}.*scenWeights{phaseIdx}(scenIdx);
                mSqrtBetaDoseExp{phaseIdx} = mSqrtBetaDoseExp{phaseIdx} + scenarioSqrtBetaDose{1}.*scenWeights{phaseIdx}(scenIdx);

            end

            lineLengthScen = lineLengthScen + fprintf('done.\n');
            lineLengthScen = lineLengthScen + fprintf('\t\t\t\tAccumulating omega for struct: ');
            
            structCounter = 0;
            
            for structIdx=structsToInclude
                structCounter = structCounter+1;
                lineLengthStruct = fprintf('%d/%d',structCounter,numel(structsToInclude));

                % Get all voxels in current structure on all ct scenarios
                currStructVoxels = [];
                for tmpScenIdx=1:nScenariosInPhase
                   currMeta = scenariosMetaInPhase(tmpScenIdx);
                   currStructVoxels = [currStructVoxels; cst{structIdx,4}{currMeta.ctScenIdx}];
                end
                currStructVoxels = unique(currStructVoxels);

                currStructOmega = currPhaseOmega{structIdx};
                currDist = gpuArray(scenarioDistribution{1}(currStructVoxels,:));

                currStructOmega = currStructOmega + gather(currDist'*currDist).*scenWeights{phaseIdx}(scenIdx);
                wait(gpu);

                currPhaseOmega{structIdx} = currStructOmega;

                if bioQuantities
                    % omega Alpha
                    currStructOmegaAlpha = currPhaseOmegaAlpha{structIdx};
                    currDistAlpha = gpuArray(scenarioAlphaDose{1}(currStructVoxels,:));
                    currStructOmegaAlpha = currStructOmegaAlpha + gather(currDistAlpha'*currDistAlpha).*scenWeights{phaseIdx}(scenIdx);
                    wait(gpu);
                    
                    currPhaseOmegaAlpha{structIdx} = currStructOmegaAlpha;
    
                    % Omega Beta
                    currStructOmegaBeta = currPhaseOmegaBeta{structIdx};
                    currDistSqrtBeta = gpuArray(scenarioSqrtBetaDose{1}(currStructVoxels,:));
                    currStructOmegaBeta = currStructOmegaBeta + gather(currDistSqrtBeta'*currDistSqrtBeta).*scenWeights{phaseIdx}(scenIdx);
                    wait(gpu);
                    currPhaseOmegaBeta{structIdx} = currStructOmegaBeta;

                    % EffectOmega
                    currStructOmegaEffect = currPhaseOmegaEffect{structIdx};
                    currDistBeta = gpuArray(scenarioBetaDose{1}(currStructVoxels,:));
                    currStructOmegaEffect = currStructOmegaEffect + currStructOmegaAlpha + 2 * gather(currDistAlpha' * currDistBeta) + gather(currDistBeta' * currDistBeta);
                    
                    wait(gpu);
                    currPhaseOmegaEffect{structIdx} = currStructOmegaEffect;

                end

                clear currDist;
                clear currDistAlpha;
                clear currDistSqrtBeta;
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
       
        if bioQuantities

            outputProbabilisticQuantities.mAlphaDoseOmegaExp(:,phaseIdx)         = currPhaseOmegaAlpha;
            outputProbabilisticQuantities.mSqrtBetaDoseOmegaExp(:,phaseIdx)      = currPhaseOmegaBeta;
            outputProbabilisticQuantities.mEffectOmegaHP(:,phaseIdx)             = currPhaseOmegaEffect;

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