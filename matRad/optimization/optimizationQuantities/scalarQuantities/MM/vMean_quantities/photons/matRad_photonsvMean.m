classdef matRad_photonsvMean < matRad_vMeanScalarQuantity

    properties (Constant)
        quantityName = 'photonsvMean';
        requiredSubquantities = {'photonsdOmega'};

        modality = 'photons';
        dOmegaSubQt = {'photonsdOmega'};
    end

    methods
        function this = matRad_photonsvMean(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_vMeanScalarQuantity(supArg{:});
        end
    end
end