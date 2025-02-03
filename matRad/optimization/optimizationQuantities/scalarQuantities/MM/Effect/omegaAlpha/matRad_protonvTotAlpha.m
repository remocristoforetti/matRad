classdef matRad_protonvTotAlpha < matRad_SMtotalVarianceAlpha

    properties (Constant)
        quantityName = 'protonsvTotAlpha';
        requiredSubquantities = {'protonsdOmegaAlpha'};

        modality = 'protons';
    end

    methods
        function this = matRad_protonvTotAlpha(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_SMtotalVarianceAlpha(supArg{:});
        end
    end
end