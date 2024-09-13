classdef matRad_DeltaOmega < matRad_ScalarQuantity

    properties (Constant)
        quantityName = 'dOmega';
        requiredSubquantities = {};

    end

    methods
        function this = matRad_DeltaOmega(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});
        end

        function quantityOutput = computeQuantity(~, dij,struct,w)
            quantityOutput = dij.physicalDoseOmega{struct} * w; %cellfun(@(structOmega) structOmega*w, dij. 'UniformOutput',false);
        end

        function gradientOutput = projectGradient(~,dij,~,fGrad,~)
            gradientOutput = fGrad * dij.physicalDoseOmega{struct};
        end
    end
end