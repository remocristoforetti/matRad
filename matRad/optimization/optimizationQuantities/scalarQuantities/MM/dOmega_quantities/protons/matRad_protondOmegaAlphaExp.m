classdef matRad_protondOmegaAlphaExp < matRad_FluenceScalarQuantity

    properties (Constant)
        quantityName = 'protonsdOmegaAlphaExp';
        requiredSubquantities = {};

        modality = 'protons';
        dijField = {'mAlphaDoseOmegaExp'};
    end

    methods
        function this = matRad_protondOmegaAlphaExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_FluenceScalarQuantity(supArg{:});
        end
    end
end