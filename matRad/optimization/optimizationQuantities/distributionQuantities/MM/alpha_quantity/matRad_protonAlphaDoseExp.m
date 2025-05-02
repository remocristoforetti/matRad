classdef matRad_protonAlphaDoseExp < matRad_FluenceDistributionQuantity

    properties (Constant)
        quantityName = 'protonsAlphaDoseExp';
        requiredSubquantities = {};

        modality = 'protons';
        dijField = {'mAlphaDoseExp'};
    end

    methods
        function this = matRad_protonAlphaDoseExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_FluenceDistributionQuantity(supArg{:});
        end
    end
end