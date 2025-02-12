classdef matRad_protonVoxelVarianceExp < matRad_VarianceQuantity

    properties (Constant)
        quantityName = 'protonsVoxelVarianceExp';
        requiredSubquantities = {'protonsvMeanExp', 'protonsMeanDoseExp'};

        modality = 'protons';
        vMeanQuantity = 'protonsvMeanExp';
        meanQuantity = 'protonsMeanDoseExp';
    end

    methods
        function this = matRad_protonVoxelVarianceExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_VarianceQuantity(supArg{:});
        end
    end
end