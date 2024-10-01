classdef matRad_DeltaOmegaBeta < matRad_ScalarQuantity

    properties (Constant)
        quantityName = 'dOmegaBeta';
        requiredSubquantities = {};

    end

    methods
        function this = matRad_DeltaOmegaBeta(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});
        end

        function quantityOutput = computeQuantity(~, dij,struct,w)
            quantityOutput = dij.mSqrtBetaDoseOmega{struct} * w; %cellfun(@(structOmega) structOmega*w, dij. 'UniformOutput',false);
        end

        function gradientOutput = projectGradient(~,dij,~,fGrad,~)
            gradientOutput = fGrad{struct} * dij.mSqrtBetaDoseOmega{struct};
        end
    end
end