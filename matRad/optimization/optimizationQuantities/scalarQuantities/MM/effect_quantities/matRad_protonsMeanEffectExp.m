classdef matRad_protonsMeanEffectExp < matRad_MeanEffect

    properties (Constant)
        quantityName = 'protonsMeanEffectExp';
        requiredSubquantities = {'protonsMeanAlphaDoseExp', 'protonsvMeanSqrtBetaExp'};

        alphaSubQuantity = 'protonsMeanAlphaDoseExp';
        betaSubQuantity = 'protonsvMeanSqrtBetaExp';
    end

    methods
        function this = matRad_protonsMeanEffectExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_MeanEffect(supArg{:});
        end
    end
end