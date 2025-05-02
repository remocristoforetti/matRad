classdef matRad_photonsdOmegaEffectHP < matRad_FluenceScalarQuantity

    properties (Constant)
        quantityName = 'photonsdOmegaEffectHP';
        requiredSubquantities = {};

        modality = 'photons';
        dijField = {'mEffectOmegaHP'};
    end

    methods
        function this = matRad_photonsdOmegaEffectHP(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_FluenceScalarQuantity(supArg{:});
        end
    end
end