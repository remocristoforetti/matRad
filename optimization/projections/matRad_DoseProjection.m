classdef matRad_DoseProjection < matRad_BackProjection
% matRad_DoseProjection class to compute physical dose during optimization
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
        function obj = matRad_DoseProjection()
            
        end
    end
    
    methods 
        function d = computeSingleScenario(~,dij,scen,w)
            if ~isempty(dij.physicalDose{scen})
                d = dij.physicalDose{scen}*w;
            else
                d = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end 
        end
        
        function [dExp,dOmegaV, vTot] = computeSingleScenarioProb(~,dij,scen,w)
            if ~isempty(dij.physicalDoseExp{scen})
                dExp = dij.physicalDoseExp{scen}*w;
                
                selectedStructs = find(~cellfun(@isempty, dij.physicalDoseOmega(:,scen)));

                dOmegaV = cell(size(dij.physicalDoseOmega,1),1);
                vTot = cell(size(dij.physicalDoseOmega,1),1);

                
                dOmegaV(selectedStructs) = arrayfun(@(i) dij.physicalDoseOmega{i,scen}*w, selectedStructs, 'UniformOutput',false);
                vTot(selectedStructs)= arrayfun(@(i) w'*dOmegaV{i}, selectedStructs, 'UniformOutput',false);
                
                % tic
                % for i = 1:size(dij.physicalDoseOmega,1)
                %     if ~isempty(dij.physicalDoseOmega{i,scen})
                %         dOmegaV{i,1} = dij.physicalDoseOmega{i,scen} * w;
                %         vTot{i,1} = w'*dOmegaV{i,1};
                %     else
                %         dOmegaV{i,1} = [];
                %         vTot{i,1}    = [];
                %     end
                % end 
                % toc
            else
                dExp = [];
                dOmegaV = [];
                vTot    = [];
            end             
        end
        
        function wGrad = projectSingleScenarioGradient(~,dij,doseGrad,scen,~)

            if ~isempty(dij.physicalDose{scen})
                wGrad = (doseGrad{scen}' * dij.physicalDose{scen})';
            else
                wGrad = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end
        end
        
        function wGrad = projectSingleScenarioGradientProb(~,dij,dExpGrad,dOmegaVgrad,scen,~)
            if ~isempty(dij.physicalDoseExp{scen})
                wGrad = (dExpGrad{scen}' * dij.physicalDoseExp{scen})';
                wGrad = wGrad + 2 * dOmegaVgrad{scen};

            else
                wGrad = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end
        end
    end
end

