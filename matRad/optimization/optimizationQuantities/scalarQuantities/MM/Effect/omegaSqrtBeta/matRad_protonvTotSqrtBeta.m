classdef matRad_protonvTotSqrtBeta < matRad_SMtotalVarianceSqrtBeta

    properties (Constant)
        quantityName = 'protonsvTotSqrtBeta';
        requiredSubquantities = {'protonsdOmegaSqrtBeta'};

        modality = 'protons';
    end

    methods
        function this = matRad_protonvTotSqrtBeta(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_SMtotalVarianceSqrtBeta(supArg{:});
        end
    end
end