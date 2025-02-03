classdef matRad_photonvTotSqrtBeta < matRad_SMtotalVarianceSqrtBeta

    properties (Constant)
        quantityName = 'photonsvTotSqrtBeta';
        requiredSubquantities = {'photonsdOmegaSqrtBeta'};

        modality = 'photons';
    end

    methods
        function this = matRad_photonvTotSqrtBeta(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_SMtotalVarianceSqrtBeta(supArg{:});
        end
    end
end