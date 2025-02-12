classdef matRad_protondOmegaAlpha < matRad_FluenceScalarQuantity

    properties (Constant)
        quantityName = 'protonsdOmegaAlpha';
        requiredSubquantities = {};

        modality = 'protons';
        dijField = {'mAlphaDoseOmega'};
    end

    methods
        function this = matRad_protondOmegaAlpha(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_FluenceScalarQuantity(supArg{:});
        end
    end
end