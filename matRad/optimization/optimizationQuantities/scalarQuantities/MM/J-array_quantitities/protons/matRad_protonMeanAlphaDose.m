classdef matRad_protonMeanAlphaDose < matRad_MeanScalarQuantity

    properties (Constant)
        quantityName = 'protonsMeanAlphaDose';
        requiredSubquantities = {};

        modality = 'protons';
        dijField = {'alphaDoseJ'};
    end

    methods
        function this = matRad_protonMeanAlphaDose(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_MeanScalarQuantity(supArg{:});
        end
    end
end