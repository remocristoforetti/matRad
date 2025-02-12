classdef matRad_photonsdOmegaAlphaExp < matRad_FluenceScalarQuantity

    properties (Constant)
        quantityName = 'photonsdOmegaAlphaExp';
        requiredSubquantities = {};

        modality = 'photons';
        dijField = {'mAlphaDoseOmegaExp'};
    end

    methods
        function this = matRad_photonsdOmegaAlphaExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_FluenceScalarQuantity(supArg{:});
        end
    end
end