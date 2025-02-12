classdef matRad_protonSqrtBetaDose < matRad_FluenceDistributionQuantity

    properties (Constant)
        quantityName = 'protonsSqrtBetaDose';
        requiredSubquantities = {};

        modality = 'protons';
        dijField = {'mSqrtBetaDose'};
    end

    methods
        function this = matRad_protonSqrtBetaDose(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_FluenceDistributionQuantity(supArg{:});
        end
    end
end