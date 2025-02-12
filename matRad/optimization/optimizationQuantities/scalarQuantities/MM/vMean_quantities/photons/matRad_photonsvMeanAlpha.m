classdef matRad_photonsvMeanAlpha < matRad_vMeanScalarQuantity

    properties (Constant)
        quantityName = 'photonsvMeanAlpha';
        requiredSubquantities = {'photonsdOmegaAlpha'};

        modality = 'photons';
        dOmegaSubQt = {'photonsdOmegaAlpha'};
    end

    methods
        function this = matRad_photonsvMeanAlpha(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_vMeanScalarQuantity(supArg{:});
        end
    end
end