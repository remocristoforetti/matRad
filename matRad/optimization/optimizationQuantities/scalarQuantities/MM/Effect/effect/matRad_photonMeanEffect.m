classdef matRad_photonMeanEffect < matRad_SMmeanEffect

    properties (Constant)
        quantityName = 'photonsMeanEffect';
        requiredSubquantities = {'photonsAlphaDose', 'photonsvTotSqrtBeta'};

        modality = 'photons';
    end

    methods
        function this = matRad_photonMeanEffect(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_SMmeanEffect(supArg{:});
        end
    end
end