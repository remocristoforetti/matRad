classdef matRad_protonvMeanExp < matRad_vMeanScalarQuantity

    properties (Constant)
        quantityName = 'protonsvMeanExp';
        requiredSubquantities = {'protonsdOmegaExp'};

        modality = 'protons';
        dOmegaSubQt = {'protonsdOmegaExp'};
    end

    methods
        function this = matRad_protonvMeanExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_vMeanScalarQuantity(supArg{:});
        end
    end
end