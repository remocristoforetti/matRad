classdef matRad_ExpectedPhysicalDose < matRad_DistributionQuantity

    properties (Constant)
        quantityName = 'physicalDoseExp';
        requiredSubquantities = {};
    end

    methods
        function this = matRad_ExpectedPhysicalDose(dij)

            if nargin>0
                supArg = {dij};
            else
                supArg = {};
            end
            
            this@matRad_DistributionQuantity(supArg{:});

        end

        function quantityOutput = computeQuantity(~, dij, scen,w)
            if ~isempty(dij.physicalDoseExp{scen})
                
                quantityOutput = dij.physicalDoseExp{scen}*w;
            
            else

                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
                quantityOutput = [];
            end
        end

        function gradientOutput = projectGradient(~,dij,scen,fGrad,~)
            if ~isempty(dij.physicalDoseExp{scen})
                gradientOutput = (fGrad{scen}' * dij.physicalDoseExp{scen})';
            else
                gradientOutput = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end
        end
    end

end