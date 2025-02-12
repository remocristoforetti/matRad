classdef matRad_protonvMeanAlpha < matRad_vMeanScalarQuantity

    properties (Constant)
        quantityName = 'protonsvMeanAlpha';
        requiredSubquantities = {'protonsdOmegaAlpha'};

        modality = 'protons';
        dOmegaSubQt = {'protonsdOmegaAlpha'};
    end

    methods
        function this = matRad_protonvMeanAlpha(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_vMeanScalarQuantity(supArg{:});
        end
    end
end