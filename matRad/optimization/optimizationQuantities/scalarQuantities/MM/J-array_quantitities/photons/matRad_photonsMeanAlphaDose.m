classdef matRad_photonsMeanAlphaDose < matRad_MeanScalarQuantity

    properties (Constant)
        quantityName = 'photonsMeanAlphaDose';
        requiredSubquantities = {};

        modality = 'photons';
        dijField = {'alphaDoseJ'};
    end

    methods
        function this = matRad_photonsMeanAlphaDose(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_MeanScalarQuantity(supArg{:});
        end
    end
end