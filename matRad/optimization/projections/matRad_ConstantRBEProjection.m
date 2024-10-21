classdef matRad_ConstantRBEProjection < matRad_BackProjection
% matRad_BackProjection for optimization based on RBExDose
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    methods
        function obj = matRad_ConstantRBEProjection()
        end   
    end
    
    methods 
        function RBExD = computeSingleScenario(~,dij,scen,w)
            if ~isempty(dij.physicalDose{scen})
                RBExD = dij.physicalDose{scen} * (dij.RBE * w);
            else
                RBExD = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end 
        end
        
        function [dExp,dOmegaV,vTot] = computeSingleScenarioProb(obj,dij,scen,w)
            if ~isempty(dij.physicalDoseExp{scen})
                dExp = dij.physicalDoseExp{scen}*(dij.RBE * w);
                
                selectedStructs = obj.useStructsForOmega;

                dOmegaV = cell(size(dij.physicalDoseOmega,1),1);
                vTot = cell(size(dij.physicalDoseOmega,1),1);

                %Here omega need dij.RBE^2
                dOmegaV(selectedStructs) = arrayfun(@(i) dij.physicalDoseOmega{i,scen}*(dij.RBE * dij.RBE * w), selectedStructs, 'UniformOutput',false);
                vTot(selectedStructs)= arrayfun(@(i) w'*dOmegaV{i}, selectedStructs, 'UniformOutput',false);
                
                if any([vTot{:}]<0)
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispWarning('Negative total variance detected, this should not happen.');
                end
            else
                dExp = [];
                dOmegaV = [];
                vTot    = [];
            end             
        end
        
        function wGrad = projectSingleScenarioGradient(~,dij,doseGrad,scen,~)
            if ~isempty(dij.physicalDose{scen})
                wGrad = ((dij.RBE * doseGrad{scen})' * dij.physicalDose{scen})';
            else
                wGrad = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end
        end
        
        function wGrad = projectSingleScenarioGradientProb(~,dij,dExpGrad,dOmegaVgrad,scen,~)
            if ~isempty(dij.physicalDoseExp{scen})
                wGrad = ((dij.RBE * dExpGrad{scen})' * dij.physicalDoseExp{scen})';
                % No need here for RBE^2, it is already included in
                % dOmegaVgrad
                wGrad = wGrad + 2 * dOmegaVgrad{scen};
            else
                wGrad = [];
                matRad_cfg = MatRad_Config.instance();

                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end
        end
    end
    

end

