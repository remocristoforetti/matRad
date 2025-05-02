classdef matRad_protonsHPEffect < matRad_Effect

    properties (Constant)
        quantityName = 'protonsHPEffect';
        requiredSubquantities = {'protonsAlphaDoseExp', 'protonsSqrtBetaDoseExp'};

        alphaSubQuantity = 'protonsAlphaDoseExp';
        betaSubQuantity = 'protonsSqrtBetaDoseExp';

        modality = 'protons';
    end

    methods
        function this = matRad_protonsHPEffect(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_Effect(supArg{:});
        end
    end
end