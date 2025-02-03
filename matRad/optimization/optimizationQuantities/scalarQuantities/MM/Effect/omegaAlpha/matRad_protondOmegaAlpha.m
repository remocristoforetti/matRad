classdef matRad_protondOmegaAlpha < matRad_SMdOmegaAlpha

    properties (Constant)
        quantityName = 'protonsdOmegaAlpha';
        requiredSubquantities = {};

        modality = 'protons';
    end

    methods
        function this = matRad_protondOmegaAlpha(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_SMdOmegaAlpha(supArg{:});
        end
    end
end