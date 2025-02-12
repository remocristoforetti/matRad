classdef matRad_photonsvMeanExp < matRad_vMeanScalarQuantity

    properties (Constant)
        quantityName = 'photonsvMeanExp';
        requiredSubquantities = {'photonsdOmegaExp'};

        modality = 'photons';
        dOmegaSubQt = {'photonsdOmegaExp'};
    end

    methods
        function this = matRad_photonsvMeanExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_vMeanScalarQuantity(supArg{:});
        end
    end
end