classdef matRad_photonsVoxelVariance < matRad_VarianceQuantity

    properties (Constant)
        quantityName = 'photonsVoxelVariance';
        requiredSubquantities = {'photonsvMean', 'photonsMeanDose'};

        modality = 'photons';
        vMeanQuantity = 'photonsvMean';
        meanQuantity  = 'photonsMeanDose';
    end

    methods
        function this = matRad_photonsVoxelVariance(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_VarianceQuantity(supArg{:});
        end
    end
end