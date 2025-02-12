classdef matRad_MMmeanDoseExp < matRad_MMscalarQuantity

    properties (Constant)
        quantityName = 'MMmeanDoseExp';
        requiredSubquantities = {'protonsMeanDoseExp', 'photonsMeanDoseExp'};
        protonSubquantity = 'protonsMeanDoseExp';
        photonSubquantity = 'photonsMeanDoseExp';
    end

    methods
        function this = matRad_MMmeanDoseExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_MMscalarQuantity(supArg{:});
            
        end

    end
end