classdef matRad_photonAlphaDoseExp < matRad_FluenceDistributionQuantity

    properties (Constant)
        quantityName = 'photonsAlphaDoseExp';
        requiredSubquantities = {};

        modality = 'photons';
        dijField = {'mAlphaDoseExp'};
    end

    methods
        function this = matRad_photonAlphaDoseExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_FluenceDistributionQuantity(supArg{:});
        end
    end
end