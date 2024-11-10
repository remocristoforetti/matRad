classdef matRad_TotalVarianceAlphaSqrtBeta < matRad_ScalarQuantity

    properties (Constant)
        quantityName = 'vAlphaSqrtBeta';
        requiredSubquantities = {'dOmegaAlphaSqrtBeta'};

    end

    methods
        function this = matRad_TotalVarianceAlphaSqrtBeta(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});


        end

        function quantityOutput = computeQuantity(this, dij, struct,w)
            
            dOmegaAlphaSqrtBeta = this.subQuantities{1}.getResult(dij,w);
            quantityOutput = w' * dOmegaAlphaSqrtBeta{struct};%cellfun(@(structdOmega) w' * structdOmega, dOmega, 'UniformOutput',false);
        end

        function gradientOutput = projectGradient(this,dij,struct,fGrad,w)
%            dOmegaAlphaSqrtBeta = (dij.mAlphaSqrtBetaDoseOmega{struct} + dij.mAlphaSqrtBetaDoseOmega{struct}') * w;
            dOmegaAlphaSqrtBeta = dij.mTwiceAlphaSqrtBetaDoseOmega{struct} * w;
            gradientOutput = fGrad{struct} * dOmegaAlphaSqrtBeta;
        end
    end
end