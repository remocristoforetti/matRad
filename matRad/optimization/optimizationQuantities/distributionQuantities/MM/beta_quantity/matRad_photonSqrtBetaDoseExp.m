classdef matRad_photonSqrtBetaDoseExp < matRad_FluenceDistributionQuantity

    properties (Constant)
        quantityName = 'photonsSqrtBetaDoseExp';
        requiredSubquantities = {};

        modality = 'photons';
        dijField = {'mSqrtBetaDoseExp'};
    end

    methods
        function this = matRad_photonSqrtBetaDoseExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_FluenceDistributionQuantity(supArg{:});
        end
    end
end