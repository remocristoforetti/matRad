classdef matRad_protondOmega < matRad_FluenceScalarQuantity

    properties (Constant)
        quantityName = 'protonsdOmega';
        requiredSubquantities = {};

        modality = 'protons';
        dijField = {'physicalDoseOmega'};
    end

    methods
        function this = matRad_protondOmega(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_FluenceScalarQuantity(supArg{:});
        end
    end
end