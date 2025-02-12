classdef matRad_photonsvMeanSqrtBetaExp < matRad_vMeanScalarQuantity

    properties (Constant)
        quantityName = 'photonsvMeanSqrtBetaExp';
        requiredSubquantities = {'photonsdOmegaSqrtBetaExp'};

        modality = 'photons';
        dOmegaSubQt = {'photonsdOmegaSqrtBetaExp'};
    end

    methods
        function this = matRad_photonsvMeanSqrtBetaExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_vMeanScalarQuantity(supArg{:});
        end
    end
end