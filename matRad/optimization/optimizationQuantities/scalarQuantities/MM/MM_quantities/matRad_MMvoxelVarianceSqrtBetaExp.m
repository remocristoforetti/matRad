classdef matRad_MMvoxelVarianceSqrtBetaExp < matRad_MMscalarQuantity

    properties (Constant)
        quantityName = 'MMvoxelVarianceSqrtBetaExp';
        requiredSubquantities = {'protonsVoxelVarianceSqrtBetaExp', 'photonsVoxelVarianceSqrtBetaExp'};
        protonSubquantity = 'protonsVoxelVarianceSqrtBetaExp';
        photonSubquantity = 'photonsVoxelVarianceSqrtBetaExp';
    end

    methods
        function this = matRad_MMvoxelVarianceSqrtBetaExp(cst)
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