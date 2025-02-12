classdef matRad_photonsMeanDoseExp < matRad_MeanScalarQuantity

    properties (Constant)
        quantityName = 'photonsMeanDoseExp';
        requiredSubquantities = {};

        modality = 'photons';
        dijField = {'physicalDoseJExp'};
    end

    methods
        function this = matRad_photonsMeanDoseExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_MeanScalarQuantity(supArg{:});
        end
    end
end