classdef matRad_protonDose < matRad_FluenceDistributionQuantity

    properties (Constant)
        quantityName = 'protonsDose';
        requiredSubquantities = {};

        modality = 'protons';
        dijField = {'physicalDose'};
    end

    methods
        function this = matRad_protonDose(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_FluenceDistributionQuantity(supArg{:});
        end
    end
end