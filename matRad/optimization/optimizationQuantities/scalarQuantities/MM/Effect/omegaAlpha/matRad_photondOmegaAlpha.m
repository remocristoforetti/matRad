classdef matRad_photondOmegaAlpha < matRad_SMdOmegaAlpha

    properties (Constant)
        quantityName = 'photonsdOmegaAlpha';
        requiredSubquantities = {};

        modality = 'photons';
    end

    methods
        function this = matRad_photondOmegaAlpha(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_SMdOmegaAlpha(supArg{:});
        end
    end
end