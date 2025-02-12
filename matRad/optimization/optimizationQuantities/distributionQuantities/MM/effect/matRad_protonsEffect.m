classdef matRad_protonsEffect < matRad_Effect

    properties (Constant)
        quantityName = 'protonsEffect';
        requiredSubquantities = {'protonsAlphaDose', 'protonsSqrtBetaDose'};

        alphaSubQuantity = 'protonsAlphaDose';
        betaSubQuantity = 'protonsSqrtBetaDose';
    end

    methods
        function this = matRad_protonsEffect(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_Effect(supArg{:});
        end
    end
end