classdef matRad_SqrtBetaDose < matRad_OptimizationQuantity

    properties (Constant)
        quantityName = 'SqrtBetaDose';
        requiredSubquantities = {};
    end

    properties
        parameter;
    end
    methods
        function this = matRad_SqrtBetaDose()
            this@matRad_OptimizationQuantity();
        end

        function  quantityOutput = computeQuantity(~, dij,scen,w)
            quantityOutput = dij.mSqrtBetaDose{scen}*w;
        
        end


        function  gradientOutput = projectGradient(this, dij,scen,fGrad,~)
            %if exist('fGrad', 'var') && ~isempty(fGrad)
            %    gradientOutput = fGrad{scen}' * dij.mSqrtBetaDose{scen};
            %else
                gradientOutput = dij.mSqrtBetaDose{scen};

            %end
        end
    end
end