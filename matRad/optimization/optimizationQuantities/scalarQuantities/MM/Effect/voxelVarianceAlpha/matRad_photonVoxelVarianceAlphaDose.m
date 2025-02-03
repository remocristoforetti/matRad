classdef matRad_photonVoxelVarianceAlphaDose < matRad_SMvoxelVarianceAlphaDose

    properties (Constant)
        quantityName = 'photonsVoxelVarianceAlphaDose';
        requiredSubquantities = {'photonsAlphaDose', 'photonsvTotAlpha'};

        modality = 'photons';
    end

    methods
        function this = matRad_photonVoxelVarianceAlphaDose(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_SMvoxelVarianceAlphaDose(supArg{:});
        end
    end
end