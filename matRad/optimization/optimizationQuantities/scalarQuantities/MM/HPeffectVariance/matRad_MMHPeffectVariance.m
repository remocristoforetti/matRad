classdef matRad_MMHPeffectVariance < matRad_MMscalarQuantity

    properties (Constant)
        quantityName = 'MMHPeffectVariance';
        requiredSubquantities = {'protonsHPeffectVariance', 'photonsHPeffectVariance'};
        protonSubquantity = 'protonsHPeffectVariance';
        photonSubquantity = 'photonsHPeffectVariance';
    end

    methods
        function this = matRad_MMHPeffectVariance(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_MMscalarQuantity(supArg{:});
            
        end

    end

    methods (Access = protected)
    
        function SF = setSF(~,value)
            
            SF.protons = value.protons.^2;
            SF.photons = value.photons.^2;

        end
    end
end