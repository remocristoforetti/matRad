classdef matRad_SqrtBetaDose < matRad_DistributionQuantity

    properties (Constant)
        quantityName = 'SqrtBetaDose';
        requiredSubquantities = {};
    end

    properties
        parameter;
    end
    methods
        function this = matRad_SqrtBetaDose(dij)
            if nargin>0
                supArg = {dij};
            else
                supArg = {};
            end
            
            this@matRad_DistributionQuantity(supArg{:});
        end

        function  quantityOutput = computeQuantity(~, dij,scen,w)
            quantityOutput = dij.mSqrtBetaDose{scen}*w;
        
        end


        function  gradientProjectionOutput = projectGradient(this, dij,scen,fGrad,~)
            %quantityGradient = this.getQuantityGradient(dij,w);
            gradientProjectionOutput = fGrad{scen}' * dij.mSqrtBetaDose{scen};
        end

        % function gradientOutput = computeQuantityGradient(this,dij,scen,~)
        %     % This is used for computing the quantity-specific (i-j)
        %     % gradient part, this is called by the higher level quantity if
        %     % needed to compute the gradient there
        %     gradientOutput = dij.mSqrtBetaDose{scen};
        % end
    end
end