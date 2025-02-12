classdef matRad_protonvMeanAlphaExp < matRad_vMeanScalarQuantity

    properties (Constant)
        quantityName = 'protonsvMeanAlphaExp';
        requiredSubquantities = {'protonsdOmegaAlphaExp'};

        modality = 'protons';
        dOmegaSubQt = {'protonsdOmegaAlphaExp'};
    end

    methods
        function this = matRad_protonvMeanAlphaExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_vMeanScalarQuantity(supArg{:});
        end
    end
end