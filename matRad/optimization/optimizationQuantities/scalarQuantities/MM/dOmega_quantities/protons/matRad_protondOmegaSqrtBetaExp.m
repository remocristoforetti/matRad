classdef matRad_protondOmegaSqrtBetaExp < matRad_FluenceScalarQuantity

    properties (Constant)
        quantityName = 'protonsdOmegaSqrtBetaExp';
        requiredSubquantities = {};

        modality = 'protons';
        dijField = {'mSqrtBetaDoseOmegaExp'};
    end

    methods
        function this = matRad_protondOmegaSqrtBetaExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_FluenceScalarQuantity(supArg{:});
        end
    end
end