classdef matRad_protonVoxelVarianceAlphaExp < matRad_VarianceQuantity

    properties (Constant)
        quantityName = 'protonsVoxelVarianceAlphaExp';
        requiredSubquantities = {'protonsvMeanAlphaExp', 'protonsMeanAlphaDoseExp'};

        modality = 'protons';
        vMeanQuantity = 'protonsvMeanAlphaExp';
        meanQuantity = 'protonsMeanAlphaDoseExp';
    end

    methods
        function this = matRad_protonVoxelVarianceAlphaExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_VarianceQuantity(supArg{:});
        end
    end
end