classdef matRad_photonsMeanAlphaDoseExp < matRad_MeanScalarQuantity

    properties (Constant)
        quantityName = 'photonsMeanAlphaDoseExp';
        requiredSubquantities = {};

        modality = 'photons';
        dijField = {'alphaDoseJExp'};
    end

    methods
        function this = matRad_photonsMeanAlphaDoseExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_MeanScalarQuantity(supArg{:});
        end
    end
end