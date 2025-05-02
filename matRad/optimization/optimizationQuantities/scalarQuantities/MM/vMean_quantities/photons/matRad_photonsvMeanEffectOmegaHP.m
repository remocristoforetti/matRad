classdef matRad_photonsvMeanEffectOmegaHP < matRad_vMeanScalarQuantity

    properties (Constant)
        quantityName = 'photonsvMeanEffectOmegaExp';
        requiredSubquantities = {'photonsdOmegaEffectHP'};

        modality = 'photons';
        dOmegaSubQt = {'photonsdOmegaEffectHP'};
    end

    methods
        function this = matRad_photonsvMeanEffectOmegaHP(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_vMeanScalarQuantity(supArg{:});
        end
    end
end