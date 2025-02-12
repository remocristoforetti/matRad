classdef matRad_MMmeanEffect < matRad_MMscalarQuantity

    properties (Constant)
        quantityName = 'MMmeanEffect';
        requiredSubquantities = {'protonsMeanEffect', 'photonsMeanEffect'};
        protonSubquantity = 'protonsMeanEffect';
        photonSubquantity = 'photonsMeanEffect';
    end

    methods
        function this = matRad_MMmeanEffect(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_MMscalarQuantity(supArg{:});
            
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