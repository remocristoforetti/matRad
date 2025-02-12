classdef matRad_MMvoxelVariance < matRad_MMscalarQuantity

    properties (Constant)
        quantityName = 'MMvoxelVariance';
        requiredSubquantities = {'protonsVoxelVariance', 'photonsVoxelVariance'};
        protonSubquantity = 'protonsVoxelVariance';
        photonSubquantity = 'photonsVoxelVariance';
    end

    methods
        function this = matRad_MMvoxelVariance(cst)
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