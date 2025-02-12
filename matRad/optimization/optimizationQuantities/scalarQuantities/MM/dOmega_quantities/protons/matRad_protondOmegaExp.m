classdef matRad_protondOmegaExp < matRad_FluenceScalarQuantity

    properties (Constant)
        quantityName = 'protonsdOmegaExp';
        requiredSubquantities = {};

        modality = 'protons';
        dijField = {'physicalDoseOmegaExp'};
    end

    methods
        function this = matRad_protondOmegaExp(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_FluenceScalarQuantity(supArg{:});
        end
    end
end