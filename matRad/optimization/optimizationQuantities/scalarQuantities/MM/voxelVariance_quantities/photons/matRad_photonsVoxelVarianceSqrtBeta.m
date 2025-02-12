classdef matRad_photonsVoxelVarianceSqrtBeta < matRad_VarianceQuantity

    properties (Constant)
        quantityName = 'photonsVoxelVarianceSqrtBeta';
        requiredSubquantities = {'photonsvMeanSqrtBeta', 'photonsMeanSqrtBetaDose'};

        modality = 'photons';
        vMeanQuantity = 'photonsvMeanSqrtBeta';
        meanQuantity = 'photonsMeanSqrtBetaDose';
    end

    methods
        function this = matRad_photonsVoxelVarianceSqrtBeta(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_VarianceQuantity(supArg{:});
        end
    end
end