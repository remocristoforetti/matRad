classdef matRad_protondOmegaEffectHP < matRad_FluenceScalarQuantity

    properties (Constant)
        quantityName = 'protonsdOmegaEffectHP';
        requiredSubquantities = {};

        modality = 'protons';
        dijField = {'mEffectOmegaHP'};
    end

    methods
        function this = matRad_protondOmegaEffectHP(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_FluenceScalarQuantity(supArg{:});
        end
    end
end