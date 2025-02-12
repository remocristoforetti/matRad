classdef matRad_protonvMeanSqrtBetaExp < matRad_vMeanScalarQuantity

    properties (Constant)
        quantityName = 'protonsvMeanSqrtBetaExp';
        requiredSubquantities = {'protonsdOmegaSqrtBetaExp'};

        modality = 'protons';
        dOmegaSubQt = {'protonsdOmegaSqrtBetaExp'};
    end

    methods
        function this = matRad_protonvMeanSqrtBetaExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_vMeanScalarQuantity(supArg{:});
        end
    end
end