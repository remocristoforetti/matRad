classdef matRad_photonsvMeanSqrtBeta < matRad_vMeanScalarQuantity

    properties (Constant)
        quantityName = 'photonsvMeanSqrtBeta';
        requiredSubquantities = {'photonsdOmegaSqrtBeta'};

        modality = 'photons';
        dOmegaSubQt = {'photonsdOmegaSqrtBeta'};
    end

    methods
        function this = matRad_photonsvMeanSqrtBeta(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_vMeanScalarQuantity(supArg{:});
        end
    end
end