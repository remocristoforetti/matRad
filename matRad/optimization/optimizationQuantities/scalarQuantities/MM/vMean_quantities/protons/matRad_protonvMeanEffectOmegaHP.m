classdef matRad_protonvMeanEffectOmegaHP < matRad_vMeanScalarQuantity

    properties (Constant)
        quantityName = 'protonsvMeanEffectOmegaExp';
        requiredSubquantities = {'protonsdOmegaEffectHP'};

        modality = 'protons';
        dOmegaSubQt = {'protonsdOmegaEffectHP'};
    end

    methods
        function this = matRad_protonvMeanEffectOmegaHP(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_vMeanScalarQuantity(supArg{:});
        end
    end
end