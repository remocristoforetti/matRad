classdef matRad_photonAlphaDose < matRad_SMalphaDose

    properties (Constant)
        quantityName = 'photonsAlphaDose';
        requiredSubquantities = {};

        modality = 'photons';
    end

    methods
        function this = matRad_photonAlphaDose(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_SMalphaDose(supArg{:});
        end
    end
end