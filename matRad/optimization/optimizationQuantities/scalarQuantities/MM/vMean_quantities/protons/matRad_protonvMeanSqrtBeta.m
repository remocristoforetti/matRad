classdef matRad_protonvMeanSqrtBeta < matRad_vMeanScalarQuantity

    properties (Constant)
        quantityName = 'protonsvMeanSqrtBeta';
        requiredSubquantities = {'protonsdOmegaSqrtBeta'};

        modality = 'protons';
        dOmegaSubQt = {'protonsdOmegaSqrtBeta'};
    end

    methods
        function this = matRad_protonvMeanSqrtBeta(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_vMeanScalarQuantity(supArg{:});
        end
    end
end