classdef matRad_protondOmegaSqrtBeta < matRad_SMdOmegaSqrtBeta

    properties (Constant)
        quantityName = 'protonsdOmegaSqrtBeta';
        requiredSubquantities = {};

        modality = 'protons';
    end

    methods
        function this = matRad_protondOmegaSqrtBeta(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_SMdOmegaSqrtBeta(supArg{:});
        end
    end
end