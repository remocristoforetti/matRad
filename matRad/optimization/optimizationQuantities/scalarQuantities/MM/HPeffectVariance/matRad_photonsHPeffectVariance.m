classdef matRad_photonsHPeffectVariance < matRad_HPeffectVariance

    properties (Constant)
        quantityName = 'photonsHPeffectVariance';
        requiredSubquantities = {'photonsvMeanEffectOmegaExp', 'photonsHPEffect'};

        vOmegaSubQuantity = 'photonsvMeanEffectOmegaExp';
        effectSubQuantity = 'photonsHPEffect';

        modality = 'photons';

    end

    methods

        function this = matRad_photonsHPeffectVariance(cst)

            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_HPeffectVariance(supArg{:});

        end
    end
end