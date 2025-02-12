classdef matRad_photonsdOmegaSqrtBetaExp < matRad_FluenceScalarQuantity

    properties (Constant)
        quantityName = 'photonsdOmegaSqrtBetaExp';
        requiredSubquantities = {};

        modality = 'photons';
        dijField = {'mSqrtBetaDoseOmegaExp'};
    end

    methods
        function this = matRad_photonsdOmegaSqrtBetaExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_FluenceScalarQuantity(supArg{:});
        end
    end
end