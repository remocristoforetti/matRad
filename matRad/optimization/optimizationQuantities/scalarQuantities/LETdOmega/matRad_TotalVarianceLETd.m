classdef matRad_TotalVarianceLETd < matRad_ScalarQuantity

    properties (Constant)
        quantityName = 'vLETd';
        requiredSubquantities = {'dOmegaLETd'};

    end

    methods
        function this = matRad_TotalVarianceLETd(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});


        end

        function quantityOutput = computeQuantity(this, dij, struct,w)
            
            dOmegaLETd = this.subQuantities{1}.getResult(dij,w);
            quantityOutput = w' * dOmegaLETd{struct};%cellfun(@(structdOmega) w' * structdOmega, dOmega, 'UniformOutput',false);
        end

        function gradientOutput = projectGradient(this,dij,struct,fGrad,w)
            dOmegaLETd = this.subQuantities{1}.getResult(dij,w);
            gradientOutput = 2 * fGrad{struct} * dOmegaLETd{struct}; %cellfun(@(structdOmega) fGrad * structdOmega, dOmega, 'UniformOutput', false);
        end
    end
end