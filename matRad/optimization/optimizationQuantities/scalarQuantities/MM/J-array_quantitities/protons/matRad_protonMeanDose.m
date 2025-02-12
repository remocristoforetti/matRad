classdef matRad_protonMeanDose < matRad_MeanScalarQuantity

    properties (Constant)
        quantityName = 'protonsMeanDose';
        requiredSubquantities = {};

        modality = 'protons';
        dijField = {'physicalDoseJ'};
    end

    methods
        function this = matRad_protonMeanDose(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_MeanScalarQuantity(supArg{:});
        end
    end
end