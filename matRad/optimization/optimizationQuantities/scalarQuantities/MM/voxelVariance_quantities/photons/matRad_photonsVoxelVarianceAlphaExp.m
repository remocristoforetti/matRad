classdef matRad_photonsVoxelVarianceAlphaExp < matRad_VarianceQuantity

    properties (Constant)
        quantityName = 'photonsVoxelVarianceAlphaExp';
        requiredSubquantities = {'photonsvMeanAlphaExp', 'photonsMeanAlphaDoseExp'};

        modality = 'photons';
        vMeanQuantity = 'photonsvMeanAlphaExp';
        meanQuantity = 'photonsMeanAlphaDoseExp';
    end

    methods
        function this = matRad_photonsVoxelVarianceAlphaExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_VarianceQuantity(supArg{:});
        end
    end
end