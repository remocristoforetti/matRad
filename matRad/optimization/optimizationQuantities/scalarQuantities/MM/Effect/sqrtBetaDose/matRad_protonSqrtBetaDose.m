classdef matRad_protonSqrtBetaDose < matRad_SMsqrtBetaDose

    properties (Constant)
        quantityName = 'protonsSqrtBetaDose';
        requiredSubquantities = {};

        modality = 'protons';
    end

    methods
        function this = matRad_protonSqrtBetaDose(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_SMsqrtBetaDose(supArg{:});
        end
    end
end