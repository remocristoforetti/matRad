classdef matRad_photonvTotAlpha < matRad_SMtotalVarianceAlpha

    properties (Constant)
        quantityName = 'photonsvTotAlpha';
        requiredSubquantities = {'photonsdOmegaAlpha'};

        modality = 'photons';
    end

    methods
        function this = matRad_photonvTotAlpha(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_SMtotalVarianceAlpha(supArg{:});
        end
    end
end