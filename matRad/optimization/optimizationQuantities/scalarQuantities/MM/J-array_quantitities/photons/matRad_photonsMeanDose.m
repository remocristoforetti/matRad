classdef matRad_photonsMeanDose < matRad_MeanScalarQuantity

    properties (Constant)
        quantityName = 'photonsMeanDose';
        requiredSubquantities = {};

        modality = 'photons';
        dijField = {'physicalDoseJ'};
    end

    methods
        function this = matRad_photonsMeanDose(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_MeanScalarQuantity(supArg{:});
        end
    end
end