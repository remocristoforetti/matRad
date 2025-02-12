classdef matRad_photonsdOmegaSqrtBeta < matRad_FluenceScalarQuantity

    properties (Constant)
        quantityName = 'photonsdOmegaSqrtBeta';
        requiredSubquantities = {};

        modality = 'photons';
        dijField = {'mSqrtBetaDoseOmega'};
    end

    methods
        function this = matRad_photonsdOmegaSqrtBeta(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_FluenceScalarQuantity(supArg{:});
        end
    end
end