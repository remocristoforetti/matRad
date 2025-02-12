classdef matRad_photonsVoxelVarianceExp < matRad_VarianceQuantity

    properties (Constant)
        quantityName = 'photonsVoxelVarianceExp';
        requiredSubquantities = {'photonsvMeanExp', 'photonsMeanDoseExp'};

        modality = 'photons';
        vMeanQuantity = 'photonsvMeanExp';
        meanQuantity = 'photonsMeanDoseExp';
    end

    methods
        function this = matRad_photonsVoxelVarianceExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_VarianceQuantity(supArg{:});
        end
    end
end