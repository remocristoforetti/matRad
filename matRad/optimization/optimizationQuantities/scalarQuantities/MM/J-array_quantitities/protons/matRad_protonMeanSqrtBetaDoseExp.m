classdef matRad_protonMeanSqrtBetaDoseExp < matRad_MeanScalarQuantity

    properties (Constant)
        quantityName = 'protonsMeanSqrtBetaDoseExp';
        requiredSubquantities = {};

        modality = 'protons';
        dijField = {'sqrtBetaDoseJExp'};
    end

    methods
        function this = matRad_protonMeanSqrtBetaDoseExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_MeanScalarQuantity(supArg{:});
        end
    end
end