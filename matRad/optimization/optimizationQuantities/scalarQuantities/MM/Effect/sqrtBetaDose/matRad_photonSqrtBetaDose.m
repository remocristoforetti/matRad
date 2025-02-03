classdef matRad_photonSqrtBetaDose < matRad_SMsqrtBetaDose

    properties (Constant)
        quantityName = 'photonsSqrtBetaDose';
        requiredSubquantities = {};

        modality = 'photons';
    end

    methods
        function this = matRad_photonSqrtBetaDose(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_SMsqrtBetaDose(supArg{:});
        end
    end
end