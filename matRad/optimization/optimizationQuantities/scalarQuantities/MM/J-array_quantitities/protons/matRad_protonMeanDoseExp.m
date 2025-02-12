classdef matRad_protonMeanDoseExp < matRad_MeanScalarQuantity

    properties (Constant)
        quantityName = 'protonsMeanDoseExp';
        requiredSubquantities = {};

        modality = 'protons';
        dijField = {'physicalDoseJExp'};
    end

    methods
        function this = matRad_protonMeanDoseExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_MeanScalarQuantity(supArg{:});
        end
    end
end