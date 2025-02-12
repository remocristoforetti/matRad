classdef matRad_photonsVoxelVarianceSqrtBetaExp < matRad_VarianceQuantity

    properties (Constant)
        quantityName = 'photonsVoxelVarianceSqrtBetaExp';
        requiredSubquantities = {'photonsvMeanSqrtBetaExp', 'photonsMeanSqrtBetaDoseExp'};

        modality = 'photons';
        vMeanQuantity = 'photonsvMeanSqrtBetaExp';
        meanQuantity = 'photonsMeanSqrtBetaDoseExp';
    end

    methods
        function this = matRad_photonsVoxelVarianceSqrtBetaExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_VarianceQuantity(supArg{:});
        end
    end
end