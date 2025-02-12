classdef matRad_protonsMeanEffect < matRad_MeanEffect

    properties (Constant)
        quantityName = 'protonsMeanEffect';
        requiredSubquantities = {'protonsMeanAlphaDose', 'protonsvMeanSqrtBeta'};

        alphaSubQuantity = 'protonsMeanAlphaDose';
        betaSubQuantity = 'protonsvMeanSqrtBeta';
    end

    methods
        function this = matRad_protonsMeanEffect(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_MeanEffect(supArg{:});
        end
    end
end