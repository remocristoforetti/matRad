classdef matRad_photonsMeanSqrtBetaDose < matRad_MeanScalarQuantity

    properties (Constant)
        quantityName = 'photonsMeanSqrtBetaDose';
        requiredSubquantities = {};

        modality = 'photons';
        dijField = {'sqrtBetaDoseJ'};
    end

    methods
        function this = matRad_photonsMeanSqrtBetaDose(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_MeanScalarQuantity(supArg{:});
        end
    end
end