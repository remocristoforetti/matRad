classdef matRad_TotalVarianceAlpha < matRad_ScalarQuantity

    properties (Constant)
        quantityName = 'vAlpha';
        requiredSubquantities = {'dOmegaAlpha'};

    end

    methods
        function this = matRad_TotalVarianceAlpha(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});


        end

        function quantityOutput = computeQuantity(this, dij, struct,w)
            
            dOmegaAlpha = this.subQuantities{1}.getResult(dij,w);
            quantityOutput = w' * dOmegaAlpha{struct};%cellfun(@(structdOmega) w' * structdOmega, dOmega, 'UniformOutput',false);
        end

        function gradientOutput = projectGradient(this,dij,struct,fGrad,w)
            dOmegaAlpha = this.subQuantities{1}.getResult(dij,w);
            gradientOutput = 2 * fGrad{struct} * dOmegaAlpha{struct}; %cellfun(@(structdOmega) fGrad * structdOmega, dOmega, 'UniformOutput', false);
        end
    end
end