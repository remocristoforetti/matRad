classdef matRad_protonVoxelVariance < matRad_VarianceQuantity

    properties (Constant)
        quantityName = 'protonsVoxelVariance';
        requiredSubquantities = {'protonsvMean', 'protonsMeanDose'};

        modality = 'protons';
        vMeanQuantity = 'protonsvMean';
        meanQuantity = 'protonsMeanDose';
    end

    methods
        function this = matRad_protonVoxelVariance(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_VarianceQuantity(supArg{:});
        end
    end
end