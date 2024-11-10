classdef matRad_TotalVarianceBeta < matRad_ScalarQuantity

    properties (Constant)
        quantityName = 'vBeta';
        requiredSubquantities = {'dOmegaBeta'};

    end

    methods
        function this = matRad_TotalVarianceBeta(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});


        end

        function quantityOutput = computeQuantity(this, dij, struct,w)
            
            dOmegaBeta = this.subQuantities{1}.getResult(dij,w);
            quantityOutput = w' * dOmegaBeta{struct};%cellfun(@(structdOmega) w' * structdOmega, dOmega, 'UniformOutput',false);
        end

        function gradientOutput = projectGradient(this,dij,struct,fGrad,w)
            dOmegaBeta = this.subQuantities{1}.getResult(dij,w);
            gradientOutput = 2 * fGrad{struct} * dOmegaBeta{struct}; %cellfun(@(structdOmega) fGrad * structdOmega, dOmega, 'UniformOutput', false);
        end
    end
end