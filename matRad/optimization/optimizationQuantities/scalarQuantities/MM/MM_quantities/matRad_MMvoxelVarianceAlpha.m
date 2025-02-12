classdef matRad_MMvoxelVarianceAlpha < matRad_MMscalarQuantity

    properties (Constant)
        quantityName = 'MMvoxelVarianceAlpha';
        requiredSubquantities = {'protonsVoxelVarianceAlpha', 'photonsVoxelVarianceAlpha'};
        protonSubquantity = 'protonsVoxelVarianceAlpha';
        photonSubquantity = 'photonsVoxelVarianceAlpha';
    end

    methods
        function this = matRad_MMvoxelVarianceAlpha(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_MMscalarQuantity(supArg{:});
            
        end

    end

    methods (Access = protected)
    
        function SF = setSF(this,value)
            
            SF.protons = value.protons.^2;
            SF.photons = value.photons.^2;

        end
    end
end