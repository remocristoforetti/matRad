classdef matRad_DeltaOmegaAlpha < matRad_ScalarQuantity

    properties (Constant)
        quantityName = 'dOmegaAlpha';
        requiredSubquantities = {};

    end

    methods
        function this = matRad_DeltaOmegaAlpha(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});
        end

        function quantityOutput = computeQuantity(~, dij,struct,w)
            quantityOutput = dij.mAlphaDoseOmega{struct} * w;
        end

        function gradientOutput = projectGradient(~,dij,~,fGrad,~)
            gradientOutput = fGrad{struct} * dij.mAlphaDoseOmega{struct};
        end
    end
end