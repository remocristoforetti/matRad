classdef matRad_MMHPEffect < matRad_MMdistributionQuantity

    properties (Constant)
        quantityName = 'MMHPEffect';
        requiredSubquantities = {'protonsHPEffect', 'photonsHPEffect'};
        protonSubquantity = 'protonsHPEffect';
        photonSubquantity = 'photonsHPEffect';
    end

    methods
        function this = matRad_MMHPEffect(cst)
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