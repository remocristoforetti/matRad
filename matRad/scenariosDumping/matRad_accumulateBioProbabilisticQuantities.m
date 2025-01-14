function [expDist, mAlphaDoseExp, mSqrtBetaDoseExp, mAlphaDoseOmega, mSqrtBetDoseOmega,probQuantitiesAccumulationTime] = matRad_accumulateBioProbabilisticQuantities(saveDir,ct,cst,scenariosMeta, multiScen, mode4D)
    
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

            mAlphaDoseExp    = cell(multiScen.numOfCtScen,1);
            mAlphaDoseExp(:) = {spalloc(dijTemplate.doseGrid.numOfVoxels,dijTemplate.totalNumOfBixels,1)};

            mSqrtBetaDoseExp    = cell(multiScen.numOfCtScen,1);
            mSqrtBetaDoseExp(:) = {spalloc(dijTemplate.doseGrid.numOfVoxels,dijTemplate.totalNumOfBixels,1)};

            mAlphaDoseOmega    = cell(size(cst,1),multiScen.numOfCtScen);
            mAlphaDoseOmega(:) = {spalloc(dijTemplate.totalNumOfBixels, dijTemplate.totalNumOfBixels,1)};

            mSqrtBetDoseOmega    = cell(size(cst,1),multiScen.numOfCtScen);
            mSqrtBetDoseOmega(:) = {spalloc(dijTemplate.totalNumOfBixels, dijTemplate.totalNumOfBixels,1)};


  
        case 'all'

            expDist = {spalloc(dijTemplate.doseGrid.numOfVoxels,dijTemplate.totalNumOfBixels,1)};
            
            mAlphaDoseExp    = {spalloc(dijTemplate.doseGrid.numOfVoxels,dijTemplate.totalNumOfBixels,1)};
            mSqrtBetaDoseExp = {spalloc(dijTemplate.doseGrid.numOfVoxels,dijTemplate.totalNumOfBixels,1)};

            mAlphaDoseOmega   = cell(size(cst,1),1);

            mAlphaDoseOmega(:) = {spalloc(dijTemplate.totalNumOfBixels, dijTemplate.totalNumOfBixels,1)};

            mSqrtBetDoseOmega   = cell(size(cst,1),1);
            mSqrtBetDoseOmega(:) = {spalloc(dijTemplate.totalNumOfBixels, dijTemplate.totalNumOfBixels,1)};

    end


    gpu = gpuDevice();
    
    matRad_cfg.dispInfo('Accumulating probabilistic quantities E[D] & Omega[D] ...\n');
    
    for phaseIdx=1:numel(ctAccumIx)
        lineLengthPhase = fprintf('\t\t4D-Phase %d/%d...\n',phaseIdx,numel(ctAccumIx));

        currPhaseOmegaAlpha = mAlphaDoseOmega(:,phaseIdx);
        currPhaseOmegaBeta  = mSqrtBetDoseOmega(:,phaseIdx);

        % Get scenarios meta in this psecific phase
        scenariosMetaInPhase = scenariosMeta(ctAccumIx{phaseIdx}(:,4));
        nScenariosInPhase = numel(scenariosMetaInPhase);

        % Accumulate quantities
        for scenIdx=1:nScenariosInPhase
            lineLengthScen = fprintf('\t\t\tAccumulating scenario: %d/%d\n',scenIdx,nScenariosInPhase);
            currMeta = scenariosMetaInPhase(scenIdx);
            [scenarioDistribution,scenarioAlphaDose,scenarioSqrtBetaDose] = matRad_loadScenariosFromMeta(saveDir,currMeta,dijTemplate);

            lineLengthScen = lineLengthScen + fprintf('\t\t\t\tAccumulating exp...');

            expDist{phaseIdx} = expDist{phaseIdx} + scenarioDistribution{1}.*scenWeights{phaseIdx}(scenIdx);

            mAlphaDoseExp{phaseIdx}    = mAlphaDoseExp{phaseIdx} + scenarioAlphaDose{1}.*scenWeights{phaseIdx}(scenIdx);
            mSqrtBetaDoseExp{phaseIdx} = mSqrtBetaDoseExp{phaseIdx} + scenarioSqrtBetaDose{1}.*scenWeights{phaseIdx}(scenIdx);

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
                currStructVoxels = cst{structIdx,4}{currMeta.ctScenIdx};

                % omega Alpha
                currStructOmegaAlpha = currPhaseOmegaAlpha{structIdx};
                currDistAlpha = gpuArray(scenarioAlphaDose{1}(currStructVoxels,:));

                currStructOmegaAlpha = currStructOmegaAlpha + gather(currDistAlpha'*currDistAlpha).*scenWeights{phaseIdx}(scenIdx);
                wait(gpu);
                clear currDistAlpha;

                currPhaseOmegaAlpha{structIdx} = currStructOmegaAlpha;

                % Omega Beta
                currStructOmegaBeta = currPhaseOmegaBeta{structIdx};
                currDistBeta = gpuArray(scenarioSqrtBetaDose{1}(currStructVoxels,:));

                currStructOmegaBeta = currStructOmegaBeta + gather(currDistBeta'*currDistBeta).*scenWeights{phaseIdx}(scenIdx);
                wait(gpu);
                clear currDistBeta;

                currPhaseOmegaBeta{structIdx} = currStructOmegaBeta;

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


        % for structIdx=structsToInclude
        % 
        %     % This done only for beta term
        %     currStructVoxels = cst{structIdx,4}{currMeta.ctScenIdx};
        % 
        %     currSqrtBetaDoseExp = mSqrtBetaDoseExp{phaseIdx}(currStructVoxels,:);
        %     currStructOmegaBeta = currPhaseOmegaBeta{structIdx};
        % 
        %     currPhaseOmegaBeta{structIdx} = sparse(currStructOmegaBeta - currSqrtBetaDoseExp'*currSqrtBetaDoseExp);
        % 
        %     clear currStructOmegaBeta;
        % 
        % end

        mAlphaDoseOmega(:,phaseIdx) = currPhaseOmegaAlpha;
        mSqrtBetDoseOmega(:,phaseIdx)  = currPhaseOmegaBeta;

        if phaseIdx~=numel(ctAccumIx)
            fprintf(repmat('\b',1,lineLengthPhase));
            lineLengthPhase = 0;
        end

    end
    matRad_cfg.dispInfo('done.\n');
    probQuantitiesAccumulationTime = toc;

end