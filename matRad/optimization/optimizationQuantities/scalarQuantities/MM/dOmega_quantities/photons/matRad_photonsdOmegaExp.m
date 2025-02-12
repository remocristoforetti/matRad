classdef matRad_photonsdOmegaExp < matRad_FluenceScalarQuantity

    properties (Constant)
        quantityName = 'photonsdOmegaExp';
        requiredSubquantities = {};

        modality = 'photons';
        dijField = {'physicalDoseOmegaExp'};
    end

    methods
        function this = matRad_photonsdOmegaExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_FluenceScalarQuantity(supArg{:});
        end
    end
end