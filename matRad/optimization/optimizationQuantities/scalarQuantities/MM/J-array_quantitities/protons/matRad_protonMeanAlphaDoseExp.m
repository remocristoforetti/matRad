classdef matRad_protonMeanAlphaDoseExp < matRad_MeanScalarQuantity

    properties (Constant)
        quantityName = 'protonsMeanAlphaDoseExp';
        requiredSubquantities = {};

        modality = 'protons';
        dijField = {'alphaDoseJExp'};
    end

    methods
        function this = matRad_protonMeanAlphaDoseExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_MeanScalarQuantity(supArg{:});
        end
    end
end