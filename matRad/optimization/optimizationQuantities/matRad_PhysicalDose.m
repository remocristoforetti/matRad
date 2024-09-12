classdef matRad_PhysicalDose < matRad_OptimizationQuantity

    properties (Constant)
        quantityName = 'physicalDose';
        requiredSubquantities = {};
    end

    methods
        function this = matRad_PhysicalDose()
            this@matRad_OptimizationQuantity();
        end

        function quantityOutput = computeQuantity(~, dij, scen,w)
            if ~isempty(dij.physicalDose{scen})
                
                quantityOutput = dij.physicalDose{scen}*w;
            
            else

                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
                quantityOutput = [];
            end
        end

        function gradientOutput = projectGradient(~,dij,scen,fGrad,~)
            if ~isempty(dij.physicalDose{scen})
                gradientOutput = (fGrad{scen}' * dij.physicalDose{scen})';
            else
                gradientOutput = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end
        end
    end

end