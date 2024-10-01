function weightGradient = matRad_objectiveGradient(optiProb,w,dij,cst)
% matRad IPOPT callback: gradient function for inverse planning
% supporting mean dose objectives, EUD objectives, squared overdosage,
% squared underdosage, squared deviation and DVH objectives
%
% call
%   g = matRad_gradFuncWrapper(optiProb,w,dij,cst)
%
% input
%   optiProb: option struct defining the type of optimization
%   w:       bixel weight vector
%   dij:     dose influence matrix
%   cst:     matRad cst struct
%
% output
%   g: gradient of objective function
%
% References
%   [1] http://www.sciencedirect.com/science/article/pii/S0958394701000577
%   [2] http://www.sciencedirect.com/science/article/pii/S0360301601025858
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 the matRad development team.
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

% get current dose / effect / RBExDose vector
optiProb.BP.compute(dij,w);
d = optiProb.BP.GetResult();

% also get probabilistic quantities (nearly no overhead if empty)
%[dExp,dOmega, vTot] = optiProb.BP.GetResultProb();

% get the used scenarios
useScen  = optiProb.BP.scenarios;
scenProb = optiProb.BP.scenarioProb;
useNominalCtScen = optiProb.BP.nominalCtScenarios;

% retrieve matching 4D scenarios
fullScen      = cell(ndims(d),1);
[fullScen{:}] = ind2sub(size(d),useScen);
% fullScen      = cell(ndims(dExp),1);
% [fullScen{:}] = ind2sub(size(dExp),useScen);

contourScen   = fullScen{1};

doseGradient          = cell(size(dij.physicalDose));
doseGradient(useScen) = {zeros(dij.doseGrid.numOfVoxels,1)};

%For probabilistic optimization
[dExp,dOmega,vTot] = optiProb.BP.GetResultProb();



if ~isempty(dExp)
    nonEmptyExp = find(~cellfun(@isempty, dExp))';
else
    nonEmptyExp = [];
end

vOmega                   = cell(numel(nonEmptyExp),1); 
vOmega(nonEmptyExp) = {zeros(dij.totalNumOfBixels,1)};


doseGradientOmega = cell(numel(nonEmptyExp),1);
doseGradientOmega(nonEmptyExp) = {zeros(dij.doseGrid.numOfVoxels,1)};

%For COWC
f_COWC = zeros(size(dij.physicalDose));

% compute objective function for every VOI.
for  i = 1:size(cst,1)
   
    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )
        
        % loop over the number of constraints and objectives for the current VOI
        for j = 1:numel(cst{i,6})
            
            %Get current optimization function
            objective = cst{i,6}{j};
            
            % only perform gradient computations for objectives
            if isa(objective,'DoseObjectives.matRad_DoseObjective')
                
                % retrieve the robustness type
                robustness = objective.robustness;
                
                % rescale dose parameters to biological optimization quantity if required
                objective = optiProb.BP.setBiologicalDosePrescriptions(objective,cst{i,5}.alphaX,cst{i,5}.betaX);
                
                if objective.isActive
                    switch robustness
                        case 'none' % if conventional opt: just sum objectiveectives of nominal dose
                            for s = useNominalCtScen
                                ixScen = s;    %useScen(s);
                                ixContour = s; %contourScen(s);
                                d_i = d{ixScen}(cst{i,4}{ixContour});
                                %add to dose gradient
                                doseGradient{ixScen}(cst{i,4}{ixContour}) = doseGradient{ixScen}(cst{i,4}{ixContour}) + objective.penalty*objective.computeDoseObjectiveGradient(d_i);
                            end
                        case 'STOCH' % perform stochastic optimization with weighted / random scenarios
                            for s = 1:numel(useScen)
                                ixScen = useScen(s);
                                ixContour = contourScen(s);
                                
                                d_i = d{ixScen}(cst{i,4}{ixContour});
                                
                                doseGradient{ixScen}(cst{i,4}{ixContour}) = doseGradient{ixScen}(cst{i,4}{ixContour}) + ...
                                    (objective.penalty*objective.computeDoseObjectiveGradient(d_i) * scenProb(s));
                                
                            end
    
                        case 'PROB' % use the expectation value and the integral variance influence matrix
                            %First check the speficic cache for probabilistic
    
                            if ~exist('doseGradientExp','var')
                                optiProb.BP.compute(dij,w);
    
                                [dExp,dOmega,vTot] = optiProb.BP.GetResultProb();
                                
                                nonEmptyExp = find(~cellfun(@isempty, dExp))';
    
                                for s=nonEmptyExp
                                    [doseGradientExp(s,1)] = {zeros(dij.doseGrid.numOfVoxels,1)};
                                end
                            end
    
    
                            nonEmptyExp = find(~cellfun(@isempty, dExp))';
    
                            if ~isequal(nonEmptyExp,useNominalCtScen)
                                totIdx = cat(1,cst{i,4}{useNominalCtScen});
                                
                                newIdx{1} = unique(totIdx);
                            else
                                newIdx = cst{i,4}(useNominalCtScen);
                            end
    
                            for s=nonEmptyExp
                                d_i = dExp{s}(newIdx{s});
                                doseGradientExp{s}(newIdx{s}) = doseGradientExp{s}(newIdx{s}) + objective.penalty*objective.computeDoseObjectiveGradient(d_i);
                            
                                %p = objective.penalty/numel(cst{i,4}{s});
                                %p = objective.penalty/numel(d_i);
    
                                %vOmega{s,1} = vOmega{s,1} + p * dOmega{i,s};
                            end
    
                        % case 'PROB' % use the expectation value and the integral variance influence matrix
                        %     %First check the speficic cache for probabilistic
                        %     if ~exist('doseGradientExp','var')
                        %         for s=useNominalCtScen
                        %             [doseGradientExp(s,1)] = {zeros(dij.doseGrid.numOfVoxels,1)};
                        %         end
                        %     end
                        %     for s=useNominalCtScen
                        %         d_i = dExp{s}(cst{i,4}{s});
                        % 
                        %         doseGradientExp{s}(cst{i,4}{s}) = doseGradientExp{s}(cst{i,4}{s}) + objective.penalty*objective.computeDoseObjectiveGradient(d_i);
                        % 
                        %         %p = objective.penalty/numel(cst{i,4}{s});
                        %         p = objective.penalty/numel(d_i);
                        % 
                        %         vOmega{s,1} = vOmega{s,1} + p * dOmega{i,s};
                        %     end
                        case 'VWWC'  % voxel-wise worst case - takes minimum dose in TARGET and maximum in OAR
                            contourIx = unique(contourScen);
                            if ~isscalar(contourIx)
                                % voxels need to be tracked through the 4D CT,
                                % not yet implemented
                                matRad_cfg.dispError('4D VWWC optimization is currently not supported');
                            end
                            
                            % prepare min/max dose vector for voxel-wise worst case
                            if ~exist('d_tmp','var')
                                d_tmp = [d{useScen}];
                            end
                            
                            d_Scen = d_tmp(cst{i,4}{contourIx},:);
                            [d_max,max_ix] = max(d_Scen,[],2);
                            [d_min,min_ix] = min(d_Scen,[],2);
                            
                            if isequal(cst{i,3},'OAR')
                                d_i = d_max;
                            elseif isequal(cst{i,3},'TARGET')
                                d_i = d_min;
                            end
                            
                            if any(isnan(d_i))
                                matRad_cfg.dispWarning('%d NaN values in gradient.',numel(isnan(d_i)));
                            end
                            
                            deltaTmp = objective.penalty*objective.computeDoseObjectiveGradient(d_i);
                            
                            for s = 1:numel(useScen)
                                ixScen = useScen(s);
                                ixContour = contourScen(s);
                                
                                if isequal(cst{i,3},'OAR')
                                    currWcIx = double(max_ix == s);
                                elseif isequal(cst{i,3},'TARGET')
                                    currWcIx = double(min_ix == s);
                                end
                                
                                doseGradient{ixScen}(cst{i,4}{ixContour}) = doseGradient{ixScen}(cst{i,4}{ixContour}) + deltaTmp.*currWcIx;
                            end
                            
                        case 'VWWC_INV'  % voxel-wise worst case - takes minimum dose in TARGET and maximum in OAR
                            contourIx = unique(contourScen);
                            if ~isscalar(contourIx)
                                % voxels need to be tracked through the 4D CT,
                                % not yet implemented
                                matRad_cfg.dispError('4D VWWC optimization is currently not supported');
                            end
                            
                            % prepare min/max dose vector for voxel-wise worst case
                            if ~exist('d_tmp','var')
                                d_tmp = [d{useScen}];
                            end
                            
                            d_Scen = d_tmp(cst{i,4}{1},:);
                            [d_max,max_ix] = max(d_Scen,[],2);
                            [d_min,min_ix] = min(d_Scen,[],2);
                            
                            if isequal(cst{i,3},'OAR')
                                d_i = d_min;
                            elseif isequal(cst{i,3},'TARGET')
                                d_i = d_max;
                            end
                            
                            if any(isnan(d_i))
                                matRad_cfg.dispWarning('%d NaN values in gradFuncWrapper.',numel(isnan(d_i)));
                            end
                            
                            deltaTmp = objective.penalty*objective.computeDoseObjectiveGradient(d_i);
                            
                            for s = 1:numel(useScen)
                                ixScen = useScen(s);
                                ixContour = contourScen(s);
                                
                                if isequal(cst{i,3},'OAR')
                                    currWcIx = double(min_ix == s);
                                elseif isequal(cst{i,3},'TARGET')
                                    currWcIx = double(max_ix == s);
                                end
                                
                                doseGradient{ixScen}(cst{i,4}{ixContour}) = doseGradient{ixScen}(cst{i,4}{ixContour}) + deltaTmp.*currWcIx;
                            end
                            
                        case 'COWC' % composite worst case consideres ovarall the worst objective function value
                            %First check the speficic cache for COWC
                            if ~exist('delta_COWC','var')
                                delta_COWC         = cell(size(doseGradient));
                                delta_COWC(useScen)    = {zeros(dij.doseGrid.numOfVoxels,1)};
                            end
                            
                            for s = 1:numel(useScen)
                                ixScen = useScen(s);
                                ixContour = contourScen(s);
                                
                                d_i = d{ixScen}(cst{i,4}{ixContour});
                                
                                f_COWC(ixScen) = f_COWC(ixScen) + objective.penalty*objective.computeDoseObjectiveFunction(d_i);
                                delta_COWC{ixScen}(cst{i,4}{ixContour}) = delta_COWC{ixScen}(cst{i,4}{ixContour}) + objective.penalty*objective.computeDoseObjectiveGradient(d_i);
                            end
                            
                        case 'OWC' % objective-wise worst case consideres the worst individual objective function value
                            %First check the speficic cache for COWC
                            f_OWC = zeros(size(doseGradient));
                            
                            if ~exist('delta_OWC','var')
                                delta_OWC = cell(size(doseGradient));
                                delta_OWC(useScen) = {zeros(dij.doseGrid.numOfVoxels,1)};
                            end
                            
                            for s = 1:numel(useScen)
                                ixScen = useScen(s);
                                ixContour = contourScen(s);
                                
                                d_i = d{ixScen}(cst{i,4}{ixContour});
                                
                                f_OWC(ixScen) = objective.penalty*objective.computeDoseObjectiveFunction(d_i);
                                
                                delta_OWC{ixScen}(cst{i,4}{ixContour}) = objective.penalty*objective.computeDoseObjectiveGradient(d_i);
                                
                            end
                              
                            switch optiProb.useMaxApprox
                                case 'logsumexp'
                                    [~,fGrad] = optiProb.logSumExp(f_OWC);
                                case 'pnorm'
                                    [~,fGrad] = optiProb.pNorm(f_OWC,numel(useScen));
                                case 'none'
                                    [~,ix] = max(f_OWC(:));
                                    fGrad = zeros(size(f_OWC));
                                    fGrad(ix) = 1;
                                case 'otherwise'
                                    matRad_cfg.dispWarning('Unknown maximum approximation desired. Using ''none'' instead.');
                                    [~,ix] = max(f_OWC(:));
                                    fGrad = zeros(size(f_OWC));
                                    fGrad(ix) = 1;
                            end
                            
                            for s = 1:numel(useScen)
                                ixScen = useScen(s);
                                ixContour = contourScen(s);
                                if fGrad(ixScen ) ~= 0
                                    doseGradient{ixScen}(cst{i,4}{ixContour}) = doseGradient{ixScen}(cst{i,4}{ixContour}) + fGrad(ixScen)*delta_OWC{ixScen}(cst{i,4}{ixContour});
                                end
     
                            end
     
                        otherwise
                            matRad_cfg.dispError('Robustness setting %s not supported!',objective.robustness);
                            
                    end  %robustness type                              
                end
            elseif isa(objective, 'OmegaObjectives.matRad_OmegaObjective')


                robustness = objective.robustness;
        
                % rescale dose parameters to biological optimization quantity if required
                %objective = optiProb.BP.setBiologicalDosePrescriptions(objective,cst{i,5}.alphaX,cst{i,5}.betaX);
                if objective.isActive
                    switch robustness
                        case 'PROB'
                            if ~exist('vTot','var') % happens if this is the first cst struct that has PROB with OmegaObjective and no DoseObjective
                                optiProb.BP.compute(dij,w);
                                [doseGradientExp(:)] = {zeros(dij.totalNumOfBixels,1)};
                                [dExp,dOmega,vTot] = optiProb.BP.GetResultProb();
                            end
                            
                            nonEmptyExp = find(~cellfun(@isempty, dExp))';
    
                            if ~isequal(nonEmptyExp,useNominalCtScen)
                                totIdx = cat(1,cst{i,4}{useNominalCtScen});
                                
                                newIdx{1} = unique(totIdx);
                            else
                                newIdx = cst{i,4}(useNominalCtScen);
                            end
                            for s=nonEmptyExp
                                %vOmega here is sum over all structures
                                d_i = dExp{s}(newIdx{s});
                                [tvGrad, tmp_doseGradientOmega] = objective.computeTotalVarianceGradient(vTot{i,s}, d_i);
                                %tvGrad = sum(tvGrad);
                                vOmega{s,1} = vOmega{s,1} + tvGrad*dOmega{i,s} .* objective.penalty;

                                doseGradientOmega{s,1}(newIdx{s}) = doseGradientOmega{s,1}(newIdx{s}) + tmp_doseGradientOmega.*objective.penalty;
                            end
                    end
                end % isActive
            end  % objective check         
        end %objective loop       

    end %empty check    
end %cst structure loop

if exist('delta_COWC','var')   
    switch optiProb.useMaxApprox
        case 'logsumexp'
            [~,fGrad] = optiProb.logSumExp(f_COWC);
        case 'pnorm'
            [~,fGrad] = optiProb.pNorm(f_COWC,numel(useScen));
        case 'none'
            [~,ixCurrWC] = max(f_COWC(:));
            fGrad = zeros(size(f_COWC));
            fGrad(ixCurrWC) = 1;
        case 'otherwise'
            matRad_cfg.dispWarning('Unknown maximum approximation desired. Using ''none'' instead.');
            [~,ixCurrWC] = max(f_COWC(:));
            fGrad = zeros(size(f_COWC));
            fGrad(ixCurrWC) = 1;
    end
    
    for s = 1:numel(useScen)
        ixScen = useScen(s);
        if fGrad(ixScen) ~= 0
            doseGradient{ixScen} = doseGradient{ixScen} + fGrad(ixScen)*delta_COWC{ixScen};
        end
    end
end

weightGradient = zeros(dij.totalNumOfBixels,1);

optiProb.BP.computeGradient(dij,doseGradient,w);
g = optiProb.BP.GetGradient();

for s = 1:numel(useScen)
  weightGradient = weightGradient + g{useScen(s)};
end

%nonEmptyOmega = ~(cellfun(@isempty, vOmega));
%nonZerosOmega = arrayfun(@(scen) nnz(vOmega{scen}),nonEmptyOmega);
%if any(nonEmptyOmega) && all(nonZerosOmega>0)
if exist('doseGradientExp', 'var')
    if all(cellfun(@(x) sum(x) == 0, doseGradientOmega))
        doseGradientOmega = cell(numel(nonEmptyExp),1);
    end
    optiProb.BP.computeGradientProb(dij,doseGradientExp,doseGradientOmega,vOmega,w);
    gProb = optiProb.BP.GetGradientProb();
    
    %Only implemented for first scenario now
    for s=nonEmptyExp
        weightGradient = weightGradient + gProb{s};
    end
end
%end


gradientChecker = 0;
if gradientChecker == 1

    f =  matRad_objectiveFunction(optiProb,w,dij,cst);
    epsilon = 1e-6;

    ix = unique(randi([dij.totalNumOfBixels],1,5));

    for i=ix

        wInit = w;
        wInit(i) = wInit(i) + epsilon;
        fDel= matRad_objectiveFunction(optiProb,wInit,dij,cst);
        numGrad = (fDel - f)/epsilon;
        diff = (numGrad/weightGradient(i) - 1)*100;
        fprintf(['grad val #' num2str(i) '- rel diff numerical and analytical gradient = ' num2str(diff) '\n']);
        %fprintf([' any nan or zero for photons' num2str(sum(isnan(glog{1}))) ',' num2str(sum(~logical(glog{1}))) ' for protons: ' num2str(sum(isnan(glog{2}))) ',' num2str(sum(~logical(glog{2}))) '\n']);
    end
end