classdef matRad_MMmeanDose < matRad_MMscalarQuantity

    properties (Constant)
        quantityName = 'MMmeanDose';
        requiredSubquantities = {'protonsMeanDose', 'photonsMeanDose'};
        protonSubquantity = 'protonsMeanDose';
        photonSubquantity = 'photonsMeanDose';
    end

    methods
        function this = matRad_MMmeanDose(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_MMscalarQuantity(supArg{:});
            
        end

    end
end