classdef matRad_protonVoxelVarianceSqrtBetaDose < matRad_SMvoxelVarianceSqrtBetaDose

    properties (Constant)
        quantityName = 'protonsVoxelVarianceSqrtBetaDose';
        requiredSubquantities = {'protonsSqrtBetaDose', 'protonsvTotSqrtBeta'};

        modality = 'protons';
    end

    methods
        function this = matRad_protonVoxelVarianceSqrtBetaDose(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_SMvoxelVarianceSqrtBetaDose(supArg{:});
        end
    end
end