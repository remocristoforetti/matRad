classdef matRad_photonsMeanSqrtBetaDoseExp < matRad_MeanScalarQuantity

    properties (Constant)
        quantityName = 'photonsMeanSqrtBetaDoseExp';
        requiredSubquantities = {};

        modality = 'photons';
        dijField = {'sqrtBetaDoseJExp'};
    end

    methods
        function this = matRad_photonsMeanSqrtBetaDoseExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_MeanScalarQuantity(supArg{:});
        end
    end
end