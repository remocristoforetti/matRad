classdef matRad_photonsdOmegaAlpha < matRad_FluenceScalarQuantity

    properties (Constant)
        quantityName = 'photonsdOmegaAlpha';
        requiredSubquantities = {};

        modality = 'photons';
        dijField = {'mAlphaDoseOmega'};
    end

    methods
        function this = matRad_photonsdOmegaAlpha(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_FluenceScalarQuantity(supArg{:});
        end
    end
end