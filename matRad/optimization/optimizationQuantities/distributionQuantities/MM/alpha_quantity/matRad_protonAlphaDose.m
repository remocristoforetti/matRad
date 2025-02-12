classdef matRad_protonAlphaDose < matRad_FluenceDistributionQuantity

    properties (Constant)
        quantityName = 'protonsAlphaDose';
        requiredSubquantities = {};

        modality = 'protons';
        dijField = {'mAlphaDose'};
    end

    methods
        function this = matRad_protonAlphaDose(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_FluenceDistributionQuantity(supArg{:});
        end
    end
end