classdef matRad_photonsVoxelVarianceAlpha < matRad_VarianceQuantity

    properties (Constant)
        quantityName = 'photonsVoxelVarianceAlpha';
        requiredSubquantities = {'photonsvMeanAlpha', 'photonsMeanAlphaDose'};

        modality = 'photons';
        vMeanQuantity = 'photonsvMeanAlpha';
        meanQuantity = 'photonsMeanAlphaDose';
    end

    methods
        function this = matRad_photonsVoxelVarianceAlpha(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_VarianceQuantity(supArg{:});
        end
    end
end