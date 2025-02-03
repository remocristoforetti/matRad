classdef matRad_protonAlphaDose < matRad_SMalphaDose

    properties (Constant)
        quantityName = 'protonsAlphaDose';
        requiredSubquantities = {};

        modality = 'protons';
    end

    methods
        function this = matRad_protonAlphaDose(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_SMalphaDose(supArg{:});
        end
    end
end