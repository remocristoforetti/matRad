classdef matRad_photondOmegaSqrtBeta < matRad_SMdOmegaSqrtBeta

    properties (Constant)
        quantityName = 'photonsdOmegaSqrtBeta';
        requiredSubquantities = {};

        modality = 'photons';
    end

    methods
        function this = matRad_photondOmegaSqrtBeta(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_SMdOmegaSqrtBeta(supArg{:});
        end
    end
end