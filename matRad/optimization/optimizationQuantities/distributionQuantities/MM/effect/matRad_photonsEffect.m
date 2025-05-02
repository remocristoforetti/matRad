classdef matRad_photonsEffect < matRad_Effect

    properties (Constant)
        quantityName = 'photonsEffect';
        requiredSubquantities = {'photonsAlphaDose', 'photonsSqrtBetaDose'};

        alphaSubQuantity = 'photonsAlphaDose';
        betaSubQuantity = 'photonsSqrtBetaDose';

        modality = 'photons';
    end

    methods
        function this = matRad_photonsEffect(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_Effect(supArg{:});
        end
    end
end