classdef matRad_protonMeanSqrtBetaDose < matRad_MeanScalarQuantity

    properties (Constant)
        quantityName = 'protonsMeanSqrtBetaDose';
        requiredSubquantities = {};

        modality = 'protons';
        dijField = {'sqrtBetaDoseJ'};
    end

    methods
        function this = matRad_protonMeanSqrtBetaDose(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_MeanScalarQuantity(supArg{:});
        end
    end
end