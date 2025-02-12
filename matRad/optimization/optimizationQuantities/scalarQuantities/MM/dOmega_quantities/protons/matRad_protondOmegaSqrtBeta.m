classdef matRad_protondOmegaSqrtBeta < matRad_FluenceScalarQuantity

    properties (Constant)
        quantityName = 'protonsdOmegaSqrtBeta';
        requiredSubquantities = {};

        modality = 'protons';
        dijField = {'mSqrtBetaDoseOmega'};
    end

    methods
        function this = matRad_protondOmegaSqrtBeta(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_FluenceScalarQuantity(supArg{:});
        end
    end
end