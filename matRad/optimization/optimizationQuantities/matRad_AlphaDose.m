classdef matRad_AlphaDose < matRad_OptimizationQuantity & handle

    properties (Constant)
        quantityName = 'AlphaDose';
        requiredSubquantities = {};
    end

    properties

    end

    methods
        function this = matRad_AlphaDose()
            this@matRad_OptimizationQuantity();
        end

        function  quantityOutput = computeQuantity(~, dij, scen,w)
            quantityOutput = dij.mAlphaDose{scen}*w;
        end

        function  gradientOutput = projectGradient(~, dij, scen,fGrad,~)
            %if exist('fGrad', 'var') && ~isempty(fGrad)
            %    gradientOutput = fGrad{scen}' * dij.mAlphaDose{scen};
            %else
                gradientOutput = dij.mAlphaDose{scen};

            %end
        end
    end
end