classdef matRad_DeltaOmegaReduced < matRad_ScalarQuantity

    properties (Constant)
        quantityName = 'dOmegaReduced';
        requiredSubquantities = {};

    end

    methods
        function this = matRad_DeltaOmegaReduced(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});
        end

        function quantityOutput = computeQuantity(~, dij,struct,w)
            quantityOutput = dij.physicalDoseOmegaReduced{struct} * w; %cellfun(@(structOmega) structOmega*w, dij. 'UniformOutput',false);
        end

        function gradientOutput = projectGradient(~,dij,struct,fGrad,~)
            gradientOutput = fGrad{struct} * dij.physicalDoseOmegaReduced{struct};
        end
    end
end