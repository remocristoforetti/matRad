classdef matRad_MMEffect < matRad_MMdistributionQuantity

    properties (Constant)
        quantityName = 'MMEffect';
        requiredSubquantities = {'protonsEffect', 'photonsEffect'};
        protonSubquantity = 'protonsEffect';
        photonSubquantity = 'photonsEffect';
    end

    methods
        function this = matRad_MMEffect(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_MMdistributionQuantity(supArg{:});
            
        end
    end

    methods (Static)
        function optiFunc = setBiologicalDosePrescriptions(optiFunc,alphaX,betaX)
            doses = optiFunc.getDoseParameters();
            effect = alphaX*doses + betaX*doses.^2;
            optiFunc = optiFunc.setDoseParameters(effect);
        end
    end
end