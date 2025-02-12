classdef matRad_photonsdOmega < matRad_FluenceScalarQuantity

    properties (Constant)
        quantityName = 'photonsdOmega';
        requiredSubquantities = {};

        modality = 'photons';
        dijField = {'physicalDoseOmega'};
    end

    methods
        function this = matRad_photonsdOmega(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_FluenceScalarQuantity(supArg{:});
        end
    end
end