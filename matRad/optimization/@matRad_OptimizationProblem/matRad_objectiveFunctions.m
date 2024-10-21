function fIndv = matRad_objectiveFunctions(optiProb,w,dij,cst)
% matRad IPOPT objective function wrapper
%
% call
%   fIndv     = matRad_objectiveFunctions(optiProb,w,dij,cst)
%
% input
%   optiProb: matRad optimization problem
%   w:        beamlet/ pencil beam weight vector
%   dij:      matRad dose influence struct
%   cst:      matRad cst struct
%   scenario: index of dij scenario to consider (optional: default 1)
%
% output
%   fIndv:    m x n matrix that stores all objective function values
%             m: number of objectives
%             n: number of scenarios
%             there are some inconsistencies when using robust optimization            
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
    end    % retrieve matching 4D scenarios
    
    fullScen = cell(ndims(d),1);
    [fullScen{:}] = ind2sub(size(d),useScen);
    contourScen = fullScen{1};

    % initialize matrix with objective function values
    fIndv = [zeros(numel(optiProb.objectives),numel(useScen))];


    %individual objective functions
    for i = 1:size(optiProb.objIdx,1) %loop over objectives


        objective = optiProb.objectives{i}; 
        curObjIdx = optiProb.objIdx(i,1);

        %calculation differs based on robustness
        
        robustness = objective.robustness;
       
        if isa(objective,'DoseObjectives.matRad_DoseObjective')
            switch robustness
                case 'none' % if conventional opt: just sum objectives of nominal dose
    
                    for ixScen = useNominalCtScen
                        d_i = d{ixScen}(cst{curObjIdx,4}{useScen(1)});
                        fInd = objective.computeDoseObjectiveFunction(d_i);
                        fIndv(i,ixScen) = fInd;
                    end
                    
                case 'STOCH' % if prob opt: sum up expectation value of objectives
                    
                    for s = 1:numel(useScen)
                        ixScen = useScen(s);
                        ixContour = contourScen(s);
                        
                        d_i = d{ixScen}(cst{curObjIdx,4}{ixContour});
                        
                        fInd =  scenProb(s) * objective.computeDoseObjectiveFunction(d_i);
                        fIndv(i,s) = fInd;
                        
                    end
                    
                case 'PROB' % if prob opt: sum up expectation value of objectives TODO: CHECK FOR VALUE TO APPEND
                    if ~exist('dExp','var')
                        optiProb.BP.compute(dij,w);
                        [dExp,~,vTot] = optiProb.BP.GetResultProb();
                    end
        
                    if ~isequal(nonEmptyExp,useNominalCtScen)
                        totIdx = cat(1,cst{curObjIdx,4}{useNominalCtScen});
                        
                        newIdx{1} = unique(totIdx);
                    else
                        newIdx = cst{curObjIdx,4}(useNominalCtScen);
                    end
                    
                    f_objective = 0;
                    for s=nonEmptyExp
                        d_i = dExp{s}(newIdx{s});
                    
                        f_objective = f_objective + objective.computeDoseObjectiveFunction(d_i);
                    end
      
                   % singleObjective = [singleObjective,f_objective];
    
                    %if objective.isActive
                    fIndv(i,1) = f_objective;
                    %end
                    
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
                    
                    d_Scen = d_tmp(cst{curObjIdx,4}{contourIx},:);
                    
                    d_max = max(d_Scen,[],2);
                    d_min = min(d_Scen,[],2);
                    
                    if isequal(cst{curObjIdx,3},'OAR')
                        d_i = d_max;
                    elseif isequal(cst{curObjIdx,3},'TARGET')
                        d_i = d_min;
                    end
                    
                    fInd = objective.computeDoseObjectiveFunction(d_i);
                    fIndv(i,1) = fInd; %no different scenarios?
                    
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
                    
                    d_Scen = d_tmp(cst{curObjIdx,4}{contourIx},:);
                    d_max = max(d_Scen,[],2);
                    d_min = min(d_Scen,[],2);
                    
                    if isequal(cst{curObjIdx,3},'OAR')
                        d_i = d_min;
                    elseif isequal(cst{curObjIdx,3},'TARGET')
                        d_i = d_max;
                    end
                    
                    fInd = objective.computeDoseObjectiveFunction(d_i);
                    fIndv(i,1) = fInd;
                    
                case 'COWC'  % composite worst case consideres ovarall the worst objective function value
                    
                    for s = 1:numel(useScen)
                        ixScen = useScen(s);
                        ixContour = contourScen(s);
                        
                        d_i = d{ixScen}(cst{curObjIdx,4}{ixContour});
                        fIndv(i,s) = objective.computeDoseObjectiveFunction(d_i);
                    end
                    
                case 'OWC'   % objective-wise worst case considers the worst individual objective function value
                    
                    f_OWC = zeros(numel(useScen),1);
                    
                    for s = 1:numel(useScen)
                        ixScen    = useScen(s);
                        ixContour = contourScen(s);
                        
                        d_i = d{ixScen}(cst{curObjIdx,4}{ixContour});
                        f_OWC(s) =  objective.computeDoseObjectiveFunction(d_i);
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
                    fIndv(i,1) = fMax; %does this make sense? 
                    
                otherwise
                    matRad_cfg.dispError('Robustness setting %s not supported!',objective.robustness);
                    
            end  %robustness type


        elseif isa(objective, 'OmegaObjectives.matRad_OmegaObjective')

            switch robustness
                case 'PROB'

                    if ~exist('vTot','var') % happens if this is the first cst struct that has PROB with OmegaObjective and no DoseObjective
                        optiProb.BP.compute(dij,w);
                        [dExp,~,vTot] = optiProb.BP.GetResultProb();
                    end
                    
                    if ~isequal(nonEmptyExp,useNominalCtScen)
                        totIdx = cat(1,cst{curObjIdx,4}{useNominalCtScen});
                        
                        newIdx{1} = unique(totIdx);
                    else
                        newIdx = cst{curObjIdx,4}(useNominalCtScen);
                    end

                    f_objective = 0;
                    for s=nonEmptyExp%useNominalCtScen
                        f_objective = f_objective + objective.computeTotalVarianceObjective(vTot{curObjIdx,s}, numel(newIdx{s}));

                    end

                    %singleObjective = [singleObjective,f_objective];

                    %if objective.isActive

                        fIndv(i,1) = f_objective;

                    %end
            end
        end
    end
end
