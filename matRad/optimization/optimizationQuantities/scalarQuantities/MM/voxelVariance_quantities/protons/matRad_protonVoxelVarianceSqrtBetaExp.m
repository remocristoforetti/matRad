classdef matRad_protonVoxelVarianceSqrtBetaExp < matRad_VarianceQuantity

    properties (Constant)
        quantityName = 'protonsVoxelVarianceSqrtBetaExp';
        requiredSubquantities = {'protonsvMeanSqrtBetaExp', 'protonsMeanSqrtBetaDoseExp'};

        modality = 'protons';
        vMeanQuantity = 'protonsvMeanSqrtBetaExp';
        meanQuantity = 'protonsMeanSqrtBetaDoseExp';
    end

    methods
        function this = matRad_protonVoxelVarianceSqrtBetaExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_VarianceQuantity(supArg{:});
        end
    end
end