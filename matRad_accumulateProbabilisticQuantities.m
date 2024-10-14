function [expDist, omega,totalTime] = matRad_accumulateProbabilisticQuantities(saveDir, ct, cst,pln,mode4D, loadScenMode, scenarioCombination,wcSigma)

    matRad_cfg = MatRad_Config.instance();
    originaLogLevel = matRad_cfg.logLevel;
    testVariance = false;

    %% Input handling
    if ~exist('mode4D', 'var') || isempty(mode4D)
        mode4D = 'all';
    elseif ~any(strcmpi(mode4D,{'phase','all'}))
        matRad_cfg.dispError('4D calculation mode %s unknown, allowed is only ''phase'' or ''all''',mode4D);
    end

    if ~exist('loadScenMode', 'var') || isempty(loadScenMode)
        loadScenMode = 'all';
    elseif ischar(loadScenMode) && strcmp(loadScenMode, 'none')
        matRad_cfg.dispError('Cannot accumulate quantities if no scenarios are specified');
    
    end

    if ~exist('scenarioCombination', 'var') || isempty(scenarioCombination)
         scenarioCombination = 'rnd';
    end

     if ~exist('wcSigma', 'var') || isempty(wcSigma)
        wcSigma = 1;
     end
     
    dirFiles = dir(saveDir);

    if ischar(loadScenMode)
        switch loadScenMode
            case 'all'
                nAllScens = sum(arrayfun(@(file) contains(file.name, 'scenario'), dirFiles(3:end)));%sum(arrayfun(@(file) strcmp(file.name(1:8), 'scenario'), dirFiles(3:end)));
                scensToLoad = [1:nAllScens];
            otherwise
                matRad_cfg.dispError('Unreconized scenario mode: %s, only available is ''all'' or specific scenario numbers', loadScenMode);
        end
    else
        if isnumeric(loadScenMode)
            nAllScens = numel(loadScenMode);
            if ~isrow(loadScenMode)
                loadScenMode = loadScenMode';
            end
            scensToLoad = loadScenMode;
        end
    end

    % This is just to get the total number of scenarios that have been
    % generated and the correct indexing
    % [~,tmpTotalScen,~] = matRad_loadDijScenarios(ct,saveDir,'all', 'none');
    % totalLineaIdx = find(tmpTotalScen.scenMask);

    matRad_cfg.logLevel = 1;
    [dij, newMultScen, ~] = matRad_loadDijScenarios(ct,saveDir,loadScenMode, 'none', [], scenarioCombination, wcSigma);
    matRad_cfg.logLevel = originaLogLevel;
    % loadedScenarios = [];
    % for k=1:numel(loadedScenariosNames)
    %     tmpDotIdx = find(loadedScenariosNames{k} == '.');
    %     loadedScenarios = [loadedScenarios; str2num(loadedScenariosNames{k}(10:tmpDotIdx-1))];
    % end
    %% Accumulate quantities
    matRad_cfg.dispInfo('Calculating probabilistic quantities E[D] & Omega[D] ...\n');

    voiIx = [];
    for i = 1:size(cst,1)
        if ~isempty(cst{i,6})
            voiIx = [voiIx i];
        end
    end

    scens = find(newMultScen.scenMask);

    
    switch mode4D
        case 'phase'
            for ctIx = 1:newMultScen.numOfCtScen
                scenIx = newMultScen.linearMask(:,1) == ctIx;
                
    
                ctAccumIx{ctIx}(:,1:3) = newMultScen.linearMask(scenIx,:);
                ctAccumIx{ctIx}(:,4) = scensToLoad(scenIx);

                % get weights associated to thes scenarios;
                scenWeights{ctIx} = newMultScen.scenWeight(scenIx);
                
                %normalize weights for current ctScenario
                scenWeights{ctIx} = scenWeights{ctIx}./sum(scenWeights{ctIx});
    
                %cumScenInCt = cumScenInCt + sum(scenIx);
            end
                
        case 'all'
            ctAccumIx{1}(:,1:3) = newMultScen.linearMask(:,1:3);
            ctAccumIx{1}(:,4) = scensToLoad;
            scenWeights{1} = newMultScen.scenWeight./sum(newMultScen.scenWeight);
    end

    %Allocate quantities
    
    tic;
    if strcmpi(mode4D,'phase')
        expDist = cell(newMultScen.numOfCtScen,1);
        [expDist(1:newMultScen.numOfCtScen)] = {spalloc(dij.doseGrid.numOfVoxels,dij.totalNumOfBixels,1)};
        omega = cell(size(cst,1),newMultScen.numOfCtScen);
        
        ixTmp = cell(ndims(newMultScen.scenMask),1);
        [ixTmp{:}] = ind2sub(size(newMultScen.scenMask),scens);
        ctIxMap = ixTmp{1};        
    else
        expDist{1} = spalloc(dij.doseGrid.numOfVoxels,dij.totalNumOfBixels,1);
        omega = cell(size(cst,1),1);
        ctIxMap = ones(numel(scens),1);
    end
    [omega{voiIx,:}] = deal(zeros(dij.totalNumOfBixels));

    gpu = gpuDevice();
    
    for phaseIdx = 1:numel(ctAccumIx)
        matRad_cfg.dispInfo('\t\t4D-Phase %d/%d...\n',phaseIdx,numel(ctAccumIx));
        
        currOmega = omega(:,phaseIdx);

        %Add up Expected value
        scensInPhase = ctAccumIx{phaseIdx};
        
        nScen = size(scensInPhase,1);

        
        %currLinearScen = sub2ind(size(dij.physicalDose), scensInPhase(:,1), scensInPhase(:,2), scensInPhase(:,3));

        lineLength = 0;

        matRad_cfg.logLevel = 1;

        for sIx = 1:nScen
            % Load the scenario
            fprintf(repmat('\b',1,lineLength));
            lineLength = fprintf('\t\t\t Accumulating phase %d/%d, scenario %d/%d\n', phaseIdx,numel(ctAccumIx),sIx,nScen);

            lineLength = lineLength + fprintf('\t\t\t\t Loading distribution\n');

            [currDij,~,~] = matRad_loadDijScenarios(ct,saveDir,scensInPhase(sIx,4),[],0);

            % Accumulate expectation value
            lineLength = lineLength + fprintf('\t\t\t\t accumulating exp\n');
            expDist{phaseIdx} = expDist{phaseIdx}  + currDij.physicalDose{scensInPhase(sIx,1)} .* scenWeights{phaseIdx}(sIx);

            for v = voiIx
                ctIdxInPhase = unique(scensInPhase(:,1));
                structLineLength = fprintf('\t\t\t\t accumulating omega for struct %d/%d', v,numel(voiIx));
                % Get all voxels in current structure on all ct scenarios
                newIdxInPhase = [];
                for ctIx = ctIdxInPhase'
                    newIx = cst{v,4}{ctIx};
                    newIdxInPhase = [newIdxInPhase; newIx];
                end
                newIdxInPhase = unique(newIdxInPhase);
               
                currDist = gpuArray(currDij.physicalDose{scensInPhase(sIx,1)}(newIdxInPhase,:));

                currOmega{v} = currOmega{v} + gather(currDist'*currDist).*scenWeights{phaseIdx}(sIx);
                wait(gpu);
                clear currDist;

                if v ~= voiIx(end)      
                    fprintf(repmat('\b',1,structLineLength));
                end
            end
            lineLength = lineLength + structLineLength;
        end
        %lineLength = lineLength + fprintf('\t\t\t\t finalizing');
        matRad_cfg.logLevel = originaLogLevel;

        matRad_cfg.dispInfo('finalizing');
        % finalize the omega matrix

        for v = voiIx
            matRad_cfg.dispInfo('\t\t\t structure %d', v);
            ctIdxInPhase = unique(scensInPhase(:,1));
            newIdxInPhase = [];
            for ctIx = ctIdxInPhase'
                newIx = cst{v,4}{ctIx};
                newIdxInPhase = [newIdxInPhase; newIx];
            end
            newIdxInPhase = unique(newIdxInPhase);
%            currExpDist = gpuArray(expDist{phaseIdx}(newIdxInPhase,:));

%            currOmega{v, phaseIdx} = currOmega{v,phaseIdx} - gather(currExpDist'*currExpDist);

            currExpDist = expDist{phaseIdx}(newIdxInPhase,:);
            currOmega{v} = currOmega{v} - currExpDist'*currExpDist;

            
            %wait(gpu);
            [~,p] = chol(currOmega{v});
            if p > 0
                matRad_cfg.dispInfo(' (Matrix not positive definite)');
            end

            omega{v,phaseIdx} = sparse(currOmega{v});
            
            matRad_cfg.dispInfo('\n');
        end
        clear currOmega;
        matRad_cfg.dispInfo(' done!');


    end
    totalTime = toc;


end