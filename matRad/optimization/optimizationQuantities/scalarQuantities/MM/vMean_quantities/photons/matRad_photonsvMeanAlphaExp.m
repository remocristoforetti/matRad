classdef matRad_photonsvMeanAlphaExp < matRad_vMeanScalarQuantity

    properties (Constant)
        quantityName = 'photonsvMeanAlphaExp';
        requiredSubquantities = {'photonsdOmegaAlphaExp'};

        modality = 'photons';
        dOmegaSubQt = {'photonsdOmegaAlphaExp'};
    end

    methods
        function this = matRad_photonsvMeanAlphaExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_vMeanScalarQuantity(supArg{:});
        end
    end
end