classdef matRad_DeltaOmegaLETd < matRad_ScalarQuantity

    properties (Constant)
        quantityName = 'dOmegaLETd';
        requiredSubquantities = {};

    end

    methods
        function this = matRad_DeltaOmegaLETd(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});
        end

        function quantityOutput = computeQuantity(~, dij,struct,w)
            quantityOutput = dij.mLETdOmega{struct} * w; %cellfun(@(structOmega) structOmega*w, dij. 'UniformOutput',false);
        end

        function gradientOutput = projectGradient(~,dij,~,fGrad,~)
            gradientOutput = fGrad{struct} * dij.mLETdOmega{struct};
        end
    end
end