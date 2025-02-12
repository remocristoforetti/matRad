classdef matRad_MMmeanEffectExp < matRad_MMscalarQuantity

    properties (Constant)
        quantityName = 'MMmeanEffectExp';
        requiredSubquantities = {'protonsMeanEffectExp', 'photonsMeanEffectExp'};
        protonSubquantity = 'protonsMeanEffectExp';
        photonSubquantity = 'photonsMeanEffectExp';
    end

    methods
        function this = matRad_MMmeanEffectExp(cst)
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