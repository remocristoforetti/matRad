classdef matRad_MMvoxelVarianceSqrtBeta < matRad_MMscalarQuantity

    properties (Constant)
        quantityName = 'MMvoxelVarianceSqrtBeta';
        requiredSubquantities = {'protonsVoxelVarianceSqrtBeta', 'photonsVoxelVarianceSqrtBeta'};
        protonSubquantity = 'protonsVoxelVarianceSqrtBeta';
        photonSubquantity = 'photonsVoxelVarianceSqrtBeta';
    end

    methods
        function this = matRad_MMvoxelVarianceSqrtBeta(cst)
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