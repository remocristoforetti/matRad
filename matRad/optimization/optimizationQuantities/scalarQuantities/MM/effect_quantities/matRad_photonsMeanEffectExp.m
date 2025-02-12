classdef matRad_photonsMeanEffectExp < matRad_MeanEffect

    properties (Constant)
        quantityName = 'photonsMeanEffectExp';
        requiredSubquantities = {'photonsMeanAlphaDoseExp', 'photonsvMeanSqrtBetaExp'};

        alphaSubQuantity = 'photonsMeanAlphaDoseExp';
        betaSubQuantity = 'photonsvMeanSqrtBetaExp';
    end

    methods
        function this = matRad_photonsMeanEffectExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_MeanEffect(supArg{:});
        end
    end
end