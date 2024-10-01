classdef matRad_DeltaOmegaAlphaSqrtBeta < matRad_ScalarQuantity

    properties (Constant)
        quantityName = 'dOmegaAlphaSqrtBeta';
        requiredSubquantities = {};

    end

    methods
        function this = matRad_DeltaOmegaAlphaSqrtBeta(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});
        end

        function quantityOutput = computeQuantity(~, dij,struct,w)
            quantityOutput = dij.mAlphaSqrtBetaDoseOmega{struct} * w; %cellfun(@(structOmega) structOmega*w, dij. 'UniformOutput',false);
        end

        function gradientOutput = projectGradient(~,dij,~,fGrad,~)
            gradientOutput = fGrad{struct} * dij.mAlphaSqrtBetaDoseOmega{struct};
        end
    end
end