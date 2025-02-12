classdef matRad_protonVoxelVarianceSqrtBeta < matRad_VarianceQuantity

    properties (Constant)
        quantityName = 'protonsVoxelVarianceSqrtBeta';
        requiredSubquantities = {'protonsvMeanSqrtBeta', 'protonsMeanSqrtBetaDose'};

        modality = 'protons';
        vMeanQuantity = 'protonsvMeanSqrtBeta';
        meanQuantity = 'protonsMeanSqrtBetaDose';
    end

    methods
        function this = matRad_protonVoxelVarianceSqrtBeta(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_VarianceQuantity(supArg{:});
        end
    end
end