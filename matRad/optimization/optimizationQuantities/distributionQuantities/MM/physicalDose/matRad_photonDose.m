classdef matRad_photonDose < matRad_FluenceDistributionQuantity

    properties (Constant)
        quantityName = 'photonsDose';
        requiredSubquantities = {};

        modality = 'photons';
        dijField = {'physicalDose'};
    end

    methods
        function this = matRad_photonDose(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_FluenceDistributionQuantity(supArg{:});
        end

    end
end