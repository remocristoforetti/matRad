classdef matRad_protonMeanEffect < matRad_SMmeanEffect

    properties (Constant)
        quantityName = 'protonsMeanEffect';
        requiredSubquantities = {'protonsAlphaDose', 'protonsvTotSqrtBeta'};

        modality = 'protons';
    end

    methods
        function this = matRad_protonMeanEffect(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_SMmeanEffect(supArg{:});
        end
    end
end