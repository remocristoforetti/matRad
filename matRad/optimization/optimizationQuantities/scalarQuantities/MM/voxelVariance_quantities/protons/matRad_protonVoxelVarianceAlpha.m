classdef matRad_protonVoxelVarianceAlpha < matRad_VarianceQuantity

    properties (Constant)
        quantityName = 'protonsVoxelVarianceAlpha';
        requiredSubquantities = {'protonsvMeanAlpha', 'protonsMeanAlphaDose'};

        modality = 'protons';
        vMeanQuantity = 'protonsvMeanAlpha';
        meanQuantity = 'protonsMeanAlphaDose';
    end

    methods
        function this = matRad_protonVoxelVarianceAlpha(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_VarianceQuantity(supArg{:});
        end
    end
end