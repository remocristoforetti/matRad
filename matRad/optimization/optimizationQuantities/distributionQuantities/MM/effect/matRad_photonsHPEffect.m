classdef matRad_photonsHPEffect < matRad_Effect

    properties (Constant)
        quantityName = 'photonsHPEffect';
        requiredSubquantities = {'photonsAlphaDoseExp', 'photonsSqrtBetaDoseExp'};

        alphaSubQuantity = 'photonsAlphaDoseExp';
        betaSubQuantity  = 'photonsSqrtBetaDoseExp';

        modality = 'photons';
    end

    methods
        function this = matRad_photonsHPEffect(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_Effect(supArg{:});
        end
    end
end