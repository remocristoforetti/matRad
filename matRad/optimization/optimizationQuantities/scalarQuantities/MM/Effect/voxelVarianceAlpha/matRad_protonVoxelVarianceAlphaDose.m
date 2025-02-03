classdef matRad_protonVoxelVarianceAlphaDose < matRad_SMvoxelVarianceAlphaDose

    properties (Constant)
        quantityName = 'protonsVoxelVarianceAlphaDose';
        requiredSubquantities = {'protonsAlphaDose', 'protonsvTotAlpha'};

        modality = 'protons';
    end

    methods
        function this = matRad_protonVoxelVarianceAlphaDose(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_SMvoxelVarianceAlphaDose(supArg{:});
        end
    end
end