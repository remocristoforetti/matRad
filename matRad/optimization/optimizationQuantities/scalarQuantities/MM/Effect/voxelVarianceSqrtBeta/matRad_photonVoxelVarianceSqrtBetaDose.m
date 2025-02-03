classdef matRad_photonVoxelVarianceSqrtBetaDose < matRad_SMvoxelVarianceSqrtBetaDose

    properties (Constant)
        quantityName = 'photonsVoxelVarianceSqrtBetaDose';
        requiredSubquantities = {'photonsSqrtBetaDose', 'photonsvTotSqrtBeta'};

        modality = 'photons';
    end

    methods
        function this = matRad_photonVoxelVarianceSqrtBetaDose(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_SMvoxelVarianceSqrtBetaDose(supArg{:});
        end
    end
end