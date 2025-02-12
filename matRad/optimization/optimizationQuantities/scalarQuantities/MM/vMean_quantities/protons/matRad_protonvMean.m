classdef matRad_protonvMean < matRad_vMeanScalarQuantity

    properties (Constant)
        quantityName = 'protonsvMean';
        requiredSubquantities = {'protonsdOmega'};

        modality = 'protons';
        dOmegaSubQt = {'protonsdOmega'};
    end

    methods
        function this = matRad_protonvMean(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_vMeanScalarQuantity(supArg{:});
        end
    end
end