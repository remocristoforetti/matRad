function f = matRad_objectiveFunction(optiProb,w,dij,cst)
% matRad IPOPT objective function wrapper
%
% call
%   f = matRad_objectiveFuncWrapper(optiProb,w,dij,cst)
%
% input
%   optiProb: matRad optimization problem
%   w:        beamlet/ pencil beam weight vector
%   dij:      matRad dose influence struct
%   cst:      matRad cst struct
%   scenario: index of dij scenario to consider (optional: default 1)
%
% output
%   f: objective function value
%
% References
%   -
%
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

% get probabilistic quantities (nearly no overhead if empty)
[dExp,dOmega] = optiProb.BP.GetResultProb();

% get the used scenarios
useScen  = optiProb.BP.scenarios;
scenProb = optiProb.BP.scenarioProb;
useNominalCtScen = optiProb.BP.nominalCtScenarios;


if ~isempty(dExp)
    nonEmptyExp = find(~cellfun(@isempty, dExp))';
else
    nonEmptyExp =  [];
end
% retrieve matching 4D scenarios
fullScen = cell(ndims(d),1);
[fullScen{:}] = ind2sub(size(d),useScen);
% fullScen = cell(ndims(dExp),1);
% [fullScen{:}] = ind2sub(size(dExp),useScen);

contourScen = fullScen{1};

% initialize f
f = 0;

% required for COWC opt
if ~isempty(useScen)
    f_COWC = zeros(numel(useScen),1);

else
    f_COWC = 0;
end
singleObjective = [];
% compute objective function for every VOI.
for  i = 1:size(cst,1)
    
    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )
        
        % loop over the number of constraints for the current VOI
        for j = 1:numel(cst{i,6})
            
            objective = cst{i,6}{j};
            
            % only perform gradient computations for objectiveectives
            if isa(objective,'DoseObjectives.matRad_DoseObjective')
                
                % rescale dose parameters to biological optimization quantity if required
                objective = optiProb.BP.setBiologicalDosePrescriptions(objective,cst{i,5}.alphaX,cst{i,5}.betaX);

                % retrieve the robustness type
                robustness = objective.robustness;
                
                switch robustness
                    case 'none' % if conventional opt: just sum objectives of nominal dose
                        f_objective = 0;
                        for ixScen = useNominalCtScen
                            d_i = d{ixScen}(cst{i,4}{ixScen}); %{useScen(ixScen)});

                            f_objective = f_objective + objective.penalty * objective.computeDoseObjectiveFunction(d_i);

                        end
                        singleObjective = [singleObjective,f_objective];
                        
                        if objective.isActive
                            f = f + f_objective;
                        end
                    case 'STOCH' % if prob opt: sum up expectation value of objectives
                        f_objective = 0;
                        for s = 1:numel(useScen)
                            ixScen = useScen(s);
                            ixContour = contourScen(s);

                            d_i = d{ixScen}(cst{i,4}{ixContour});

                            f_objective = f_objective + scenProb(s) * objective.penalty*objective.computeDoseObjectiveFunction(d_i);
                        end
                        singleObjective = [singleObjective,f_objective];

                        if objective.isActive
                            f   = f + f_objective;
                        end

                    case 'PROB' % if prob opt: sum up expectation value of objectives
                        if ~exist('dExp','var')
                            optiProb.BP.compute(dij,w);
                            [dExp,dOmega,vTot] = optiProb.BP.GetResultProb();
                        end

                        if ~isequal(nonEmptyExp,useNominalCtScen)
                            totIdx = cat(1,cst{i,4}{useNominalCtScen});
                            
                            newIdx{1} = unique(totIdx);
                        else
                            newIdx = cst{i,4}(useNominalCtScen);
                        end
                        
                        f_objective = 0;
                        for s=nonEmptyExp
                            d_i = dExp{s}(newIdx{s});
                        
                            f_objective = f_objective + objective.penalty * objective.computeDoseObjectiveFunction(d_i);
                        end
                            singleObjective = [singleObjective,f_objective];

                            if objective.isActive
                                f = f + f_objective;
                            end
                            
                            %fprintf('struct number %d, current totF = %d, exp term = %d ', i, f, objective.penalty*objective.computeDoseObjectiveFunction(d_i));
                            %p = objective.penalty/numel(cst{i,4}{s});
                            %p = objective.penalty/numel(d_i);

                            %vOmega{s,1} = vOmega{s,1} + p * dOmega{i,s};

                        % for s=useNominalCtScen
                        %    d_i = dExp{s}(cst{i,4}{s});
                        % 
                        %    f   = f +  objective.penalty*objective.computeDoseObjectiveFunction(d_i);
                        % 
                        %    %p = objective.penalty/numel(cst{i,4}{s});
                        %    p = objective.penalty/numel(d_i);
                        % 
                        % % only one variance term per VOI
                        % %if j == 1
                        %    f = f + p * w' * dOmega{i,s};
                        % %end
                        % end
                    case 'VWWC'  % voxel-wise worst case - takes minimum dose in TARGET and maximum in OAR
                        contourIx = unique(contourScen);
                        if ~isscalar(contourIx)
                            % voxels need to be tracked through the 4D CT,
                            % not yet implemented
                            matRad_cfg.dispError('4D VWWC optimization is currently not supported');
                        end

                        % prepare min/max dose vector
                        if ~exist('d_tmp','var')
                            d_tmp = [d{useScen}];
                        end

                        d_Scen = d_tmp(cst{i,4}{contourIx},:);

                        d_max = max(d_Scen,[],2);
                        d_min = min(d_Scen,[],2);

                        if isequal(cst{i,3},'OAR')
                            d_i = d_max;
                        elseif isequal(cst{i,3},'TARGET')
                            d_i = d_min;
                        end

                        f = f + objective.penalty*objective.computeDoseObjectiveFunction(d_i);

                    case 'VWWC_INV'  %inverse voxel-wise conformitiy - consider the maximum and minimum dose in the target and optimize the dose conformity
                        contourIx = unique(contourScen);
                        if ~isscalar(contourIx)
                            % voxels need to be tracked through the 4D CT,
                            % not yet implemented
                            matRad_cfg.dispWarning('4D inverted VWWC optimization is currently not supported');
                        end

                        % prepare min/max dose vector
                        if ~exist('d_tmp','var')
                            d_tmp = [d{useScen}];
                        end

                        d_Scen = d_tmp(cst{i,4}{contourIx},:);
                        d_max = max(d_Scen,[],2);
                        d_min = min(d_Scen,[],2);

                        if isequal(cst{i,3},'OAR')
                            d_i = d_min;
                        elseif isequal(cst{i,3},'TARGET')
                            d_i = d_max;
                        end

                        f = f + objective.penalty*objective.computeDoseObjectiveFunction(d_i);

                    case 'COWC'  % composite worst case consideres ovarall the worst objective function value

                        for s = 1:numel(useScen)
                            ixScen = useScen(s);
                            ixContour = contourScen(s);

                            d_i = d{ixScen}(cst{i,4}{ixContour});

                            f_COWC(s) = f_COWC(s) + objective.penalty*objective.computeDoseObjectiveFunction(d_i);
                        end

                    case 'OWC'   % objective-wise worst case considers the worst individual objective function value

                        f_OWC = zeros(numel(useScen),1);
                        
                        for s = 1:numel(useScen)
                            ixScen    = useScen(s);
                            ixContour = contourScen(s);

                            d_i = d{ixScen}(cst{i,4}{ixContour});
                            f_OWC(s) = objective.penalty*objective.computeDoseObjectiveFunction(d_i);
                        end

                        % compute the maximum objective function value
                        switch optiProb.useMaxApprox
                            case 'logsumexp'
                                fMax = optiProb.logSumExp(f_OWC);
                            case 'pnorm'
                                fMax = optiProb.pNorm(f_OWC,numel(useScen));
                            case 'none'
                                fMax = max(f_OWC);
                            case 'otherwise'
                                matRad_cfg.dispWarning('Unknown maximum approximation desired. Using ''none'' instead.');
                                fMax = max(f_OWC);
                        end
                        f = f + fMax;
                        
                    otherwise
                        matRad_cfg.dispError('Robustness setting %s not supported!',objective.robustness);

                end  %robustness type                              
            elseif isa(objective, 'OmegaObjectives.matRad_OmegaObjective')


                %objective = optiProb.BP.setBiologicalDosePrescriptions(objective, cst{i,5}.alphaX, cst{i,5}.betaX);
                robustness = objective.robustness;

                switch robustness
                    case 'PROB'

                        if ~exist('vTot','var') % happens if this is the first cst struct that has PROB with OmegaObjective and no DoseObjective
                            optiProb.BP.compute(dij,w);
                            [dExp,~,vTot] = optiProb.BP.GetResultProb();
                        end
                        
                        if ~isequal(nonEmptyExp,useNominalCtScen)
                            totIdx = cat(1,cst{i,4}{useNominalCtScen});
                            
                            newIdx{1} = unique(totIdx);
                        else
                            newIdx = cst{i,4}(useNominalCtScen);
                        end

                        f_objective = 0;
                        for s=nonEmptyExp
                            d_i = dExp{s}(newIdx{s});
                            f_objective = f_objective + objective.penalty * objective.computeTotalVarianceObjective(vTot{i,s}, d_i);

                        end

                        singleObjective = [singleObjective,f_objective];

                        if objective.isActive

                            f = f + f_objective;

                        end
                        %fprintf(' Var =  %d\n', objective.penalty * objective.computeTotalVarianceObjective(vTot{i,s}, numel(cst{i,4}{s})));
                end
            
            
            end  % objective check         
        end %objective loop       
    end %empty check    
end %cst structure loop


%Handling the maximum of the composite worst case part
fMax = max(f_COWC(:));
if fMax > 0
    switch optiProb.useMaxApprox
        case 'logsumexp'
            fMax = optiProb.logSumExp(f_COWC);
        case 'pnorm'
            fMax = optiProb.pNorm(f_COWC,numel(useScen));
        case 'none'
            fMax = max(f_COWC);
        case 'otherwise'
            matRad_cfg.dispWarning('Unknown maximum approximation desired. Using ''none'' instead.');
            fMax = max(f_COWC);
    end
end

%if optiProb.graphicOutput.isOpen
    optiProb.graphicOutput.updateData(singleObjective, f + fMax);
    optiProb.graphicOutput.updatePlot();
%end

%Sum up max of composite worst case part

f = f + fMax;
