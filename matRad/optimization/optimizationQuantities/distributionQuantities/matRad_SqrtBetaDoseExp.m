classdef matRad_SqrtBetaDoseExp < matRad_DistributionQuantity

    properties (Constant)
        quantityName = 'SqrtBetaDoseExp';
        requiredSubquantities = {};
    end

    properties
        parameter;
    end
    methods
        function this = matRad_SqrtBetaDoseExp(dij)
            if nargin>0
                supArg = {dij};
            else
                supArg = {};
            end
            
            this@matRad_DistributionQuantity(supArg{:});
        end

        function  quantityOutput = computeQuantity(~, dij,scen,w)
            quantityOutput = dij.mSqrtBetaDoseExp{scen}*w;
        
        end


        function  gradientProjectionOutput = projectGradient(this, dij,scen,fGrad,~)
            %quantityGradient = this.getQuantityGradient(dij,w);
            gradientProjectionOutput = fGrad{scen}' * dij.mSqrtBetaDoseExp{scen};
        end

    end
end