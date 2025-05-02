classdef matRad_protonSqrtBetaDoseExp < matRad_FluenceDistributionQuantity

    properties (Constant)
        quantityName = 'protonsSqrtBetaDoseExp';
        requiredSubquantities = {};

        modality = 'protons';
        dijField = {'mSqrtBetaDoseExp'};
    end

    methods
        function this = matRad_protonSqrtBetaDoseExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_FluenceDistributionQuantity(supArg{:});
        end
    end
end