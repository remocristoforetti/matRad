classdef matRad_protonsHPeffectVariance < matRad_HPeffectVariance

    properties (Constant)
        quantityName = 'protonsHPeffectVariance';
        requiredSubquantities = {'protonsvMeanEffectOmegaExp', 'protonsHPEffect'};

        vOmegaSubQuantity = 'protonsvMeanEffectOmegaExp';
        effectSubQuantity = 'protonsHPEffect';

        modality = 'protons';
    end

    methods

        function this = matRad_protonsHPeffectVariance(cst)

            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_HPeffectVariance(supArg{:});

        end
    end
end