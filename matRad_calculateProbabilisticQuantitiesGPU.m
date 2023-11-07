function [dij] = matRad_calculateProbabilisticQuantitiesGPU(dij,cst,pln,mode4D)
% matRad helper function to compute probabilistic quantities, i.e. expected
%   dose influence and variance omega matrices for probabilistic
%   optimization
%
% call
%   [dij,cst] = matRad_calculateProbabilisticQuantities(dij,cst,pln)
%
% input
%   dij:        matRad dij struct
%   cst:        matRad cst struct (in dose grid resolution)
%   pln:        matRad pln struct
%   mode4D:     (Optional) Handling of 4D phases: 
%               - 'all'   : include 4D scen in statistic
%               - 'phase' : create statistics per phase (default)
%
% output
%   dij:        dij with added probabilistic quantities
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


matRad_cfg = MatRad_Config.instance();

testVariance = false;

if nargin < 4
    mode4D = 'phase';
elseif ~any(strcmpi(mode4D,{'phase','all'}))
    matRad_cfg.dispError('4D calculation mode %s unknown, allowed is only ''phase'' or ''all''',mode4D);
end
    

matRad_cfg.dispInfo('Calculating probabilistic quantities E[D] & Omega[D] ...\n');

if ~pln.bioParam.bioOpt
    fNames = {'physicalDose'};
else
    fNames = {'mAlphaDose','mSqrtBetaDose'};
end

% create placeholders for expected influence matrices

voiIx = [];
for i = 1:size(cst,1)
    if ~isempty(cst{i,6})
        voiIx = [voiIx i];
    end
end

scens = find(pln.multScen.scenMask);


%Set-up economic, parallelizable loops according to 4D mode

switch mode4D
    case 'phase'
        for ctIx = 1:pln.multScen.numOfCtScen
            scenIx = pln.multScen.linearMask(:,1) == ctIx;
            

            ctAccumIx{ctIx} = pln.multScen.linearMask(scenIx,:);

            % get weights associated to thes scenarios;
            scenWeights{ctIx} = pln.multScen.scenWeight(scenIx);
            
            %normalize weights for current ctScenario
            scenWeights{ctIx} = scenWeights{ctIx}./sum(scenWeights{ctIx});

            %cumScenInCt = cumScenInCt + sum(scenIx);
        end
            
    case 'all'
        ctAccumIx{1} = pln.multScen.linearMask(:,1:3);
        scenWeights{1} = pln.multScen.scenWeight./sum(pln.multScen.scenWeight);
end

%Variance-Test with vector of ones
%wTest = ones(dij.totalNumOfBixels,1);
%totVarTest = {};

        
%Create structures
for i = 1:numel(fNames)
    tic;
    matRad_cfg.dispInfo('\tE[D] & Omega[D] for ''%s'':\n',fNames{1,i});
    if strcmpi(mode4D,'phase')
        dij.([fNames{1,i} 'Exp']) = cell(pln.multScen.numOfCtScen,1);
        [dij.([fNames{1,i} 'Exp'])(1:pln.multScen.numOfCtScen)] = {spalloc(dij.doseGrid.numOfVoxels,dij.totalNumOfBixels,1)};
        dij.([fNames{1,i} 'Omega']) = cell(size(cst,1),pln.multScen.numOfCtScen);
        
        ixTmp = cell(ndims(pln.multScen.scenMask),1);
        [ixTmp{:}] = ind2sub(size(pln.multScen.scenMask),scens);
        ctIxMap = ixTmp{1};        
    else
        dij.([fNames{1,i} 'Exp']){1} = spalloc(dij.doseGrid.numOfVoxels,dij.totalNumOfBixels,1);
        dij.([fNames{1,i} 'Omega']) = cell(size(cst,1),1);
        ctIxMap = ones(numel(scens),1);
    end
    [dij.([fNames{1,i} 'Omega']){voiIx,:}] = deal(zeros(dij.totalNumOfBixels));
    
    %Now loop over the scenarios

    s = 0;

    for phaseIdx = 1:numel(ctAccumIx)
        matRad_cfg.dispInfo('\t\t4D-Phase %d/%d...\n',phaseIdx,numel(ctAccumIx));
        
        %Add up Expected value
        scensInPhase = ctAccumIx{phaseIdx};
        matRad_cfg.dispInfo('\t\tAccumulating Expected Dij ');
        
        nScen = size(scensInPhase,1);

        currLinearScen = sub2ind(size(dij.physicalDose), scensInPhase(:,1), scensInPhase(:,2), scensInPhase(:,3));

        for sIx = 1:nScen
            matRad_cfg.dispInfo('.');

            scenIx = currLinearScen(sIx);

            dij.([fNames{i} 'Exp']){phaseIdx} = dij.([fNames{i} 'Exp']){phaseIdx} + dij.([fNames{i}]){scenIx} .* scenWeights{phaseIdx}(sIx);
        end
                
        matRad_cfg.dispInfo(' done!\n');
        
        %Add up Omega
        matRad_cfg.dispInfo('\t\tAccumulating Omega Matrices:\n');
        vCnt = 0;
        
        gpu = gpuDevice();
        for v = voiIx
            
            vCnt = vCnt + 1;
            
            %remove zero rows, may give small advantage with dijs for large
            %structures
            tmpWones = ones(dij.totalNumOfBixels,1);
            tmpd = dij.([fNames{1,i} 'Exp']){phaseIdx} * tmpWones;            
            ixNz = find(tmpd > 0);
            %newIx = intersect(ixNz,cst{v,4}{});
            
            if testVariance
                testV = zeros(dij.doseGrid.numOfVoxels,size(scensInPhase,1));
            end
                        
            matRad_cfg.dispInfo('\t\t\tStructure %d/%d',vCnt,numel(voiIx));
            
            omegaCurr = dij.([fNames{1,i} 'Omega']){v,phaseIdx};
            
            ctIdxInPhase = unique(scensInPhase(:,1));
            

            % Get all voxels in current structure on all ct sceanrios
            newIdxInPhase = [];
            for ctIx = ctIdxInPhase'
                newIx = intersect(ixNz,cst{v,4}{ctIx});
                newIdxInPhase = [newIdxInPhase; newIx];
            end
            newIdxInPhase = unique(newIdxInPhase);

            for ctIx=ctIdxInPhase'

                % newIx = intersect(ixNz,cst{v,4}{ctIx});
                % newIdxInPhase = [newIdxInPhase; newIx];
                scensOnCt = find(scensInPhase(:,1)==ctIx);
                for sIx = 1:numel(scensInPhase(scensOnCt))
                    matRad_cfg.dispInfo('.');
                    
                    scenIx = currLinearScen(scensOnCt(sIx));
                    

                    dijCurr = gpuArray(dij.([fNames{1,i}]){scenIx}(newIdxInPhase,:));
    
                    if testVariance
                        testV(newIdxInPhase,scensOnCt(sIx)) = dij.([fNames{1,i}]){scenIx}(newIdxInPhase,:)*tmpWones;
                        scenWeightsForVariance(scensOnCt(sIx)) = scenWeights{phaseIdx}(scensOnCt(sIx));%pln.multScen.scenWeight(s+sIx);
                    end
                    
    
                     omegaCurr = omegaCurr + gather(dijCurr'*dijCurr).* scenWeights{phaseIdx}(scensOnCt(sIx));
    %                omegaCurr = omegaCurr + gather(dijCurr'*dijCurr).* pln.multScen.scenWeight(s+sIx);
    %                omegaCurr = omegaCurr + (dijCurr'*dijCurr).* pln.multScen.scenWeight(s+sIx);
                   
                    wait(gpu);
                    clear dijCurr omegaTmp;
                end
            
            end
            matRad_cfg.dispInfo(' finalizing...');
            
            dijCurr = dij.([fNames{1,i} 'Exp']){phaseIdx}(newIdxInPhase,:);

            omegaCurr = omegaCurr - dijCurr' * dijCurr;
            
            %
            %expPart = dij.physicalDoseExp{1}.*dij.physicalDoseExp{1};
            %varij = cumOmegaPart - expPart;
            %voxelCumulativeVarij = sum(cumOmegaPart,1) - sum(expPart,1);


            [~,p] = chol(omegaCurr);
            if p > 0
                matRad_cfg.dispInfo(' (Matrix not positive definite)');
            end

            % if p > 0
            %     matRad_cfg.dispInfo(' correcting SPD...');
            % 
            %     [Q,L] = eig(omegaCurr);
            %     L(L<0) = 0;
            %     omegaCurr = Q*L*Q';
            % 
            %     clear('Q','L')
            % end

            
            if testVariance
            	testV = var(testV,scenWeightsForVariance,2);
                testV = sum(testV);
                
                matRad_cfg.dispInfo(' Variance Test: %f (sample) v. %f (Omega)...',testV,tmpWones'*omegaCurr*tmpWones);
            end
            

            sparsity = nnz(omegaCurr)/numel(omegaCurr);
            
            %if sparsity < 0.2 
                dij.([fNames{1,i} 'Omega']){v,phaseIdx} = sparse(gather(omegaCurr));                        
            %else
            %    dij.([fNames{1,i} 'Omega']){v,ctIx} = full(gather(omegaCurr));
            %end
            
            
            clear omegaCurr;
            
            
            matRad_cfg.dispInfo(' done!\n');
        end
        
        s = s + size(scensInPhase,1);
    end

    totalTime = toc;    
    
    dij.info.OmegaCalculationTime = totalTime;
    %{
    for s = 1:numel(scens)
        matRad_cfg.dispInfo('\t\tScenario %d/%d...',s,numel(scens));
        ctIx = ctIxMap(s);
        scenIx = scens(s);
        
        
        %Add up Expected value
        dij.([fNames{1,i} 'Exp']){ctIx} = dij.([fNames{ctIx,i} 'Exp']){ctIx} + dij.([fNames{1,i}]){scenIx} .* pln.multScen.scenProb(s);
        
        
        %Add up Omega
         for v = voiIx
            dijCurr = gpuArray(dij.([fNames{1,i}]){scenIx}(cst{v,4}{ctIx},:));
            dijCurr = pln.multScen.scenProb(s) * dijCurr;
            omegaCurr = gpuArray(dij.([fNames{1,i} 'Omega']){v,ctIx});
%              dij.([fNames{1,i} 'Omega']){v,ctIx} = dij.([fNames{1,i} 'Omega']){v,ctIx} + ...
%                  (((dij.(fNames{1,i}){scenIx}(cst{v,4}{ctIx},:)' * pln.multScen.scenProb(s)) * ...
%                  (dij.(fNames{1,i}){scenIx}(cst{v,4}{ctIx},:)) * pln.multScen.scenProb(s)));
               dij.([fNames{1,i} 'Omega']){v,ctIx} = gather(omegaCurr + dijCurr'*dijCurr);
            %clear('dijCurr');
         end

         matRad_cfg.dispInfo('done!\n');
    end

    matRad_cfg.dispInfo('\tFinalizing Omega...');

    %Finalize Omega matrics
    unCtIx = unique(ctIx);    
    for ctIx = unCtIx
        for v = voiIx
            dijCurr = gpuArray(dij.([fNames{1,i} 'Exp']){ctIx}(cst{v,4}{ctIx},:));
            omegaCurr = gpuArray(dij.([fNames{1,i} 'Omega']){v,ctIx});
            %dij.([fNames{1,i} 'Omega']){v,ctIx} = dij.([fNames{1,i} 'Omega']){v,ctIx} - (dij.([fNames{1,i} 'Exp']){ctIx}(cst{v,4}{ctIx},:)' * dij.([fNames{1,i} 'Exp']){ctIx}(cst{v,4}{ctIx},:));
            dij.([fNames{1,i} 'Omega']){v,ctIx} = gather(omegaCurr - dijCurr' * dijCurr);
            %clear('dijCurr');
        end
    end  
    matRad_cfg.dispInfo('\tDone!\n');
    toc;
    
    %}
end


