classdef matRad_SMtotalVarianceSqrtBeta < matRad_ScalarQuantity

    properties (Constant)

    end

    methods
        function this = matRad_SMtotalVarianceSqrtBeta(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});
        end

        function quantityOutput = computeQuantity(this, dij, struct,w)
            
            wModality = w.(this.modality);

            dOmegaQt =  this.getSubQuantity([this.modality, 'dOmegaSqrtBeta']);
            dOmega = dOmegaQt.getResult(dij,w);
            
            quantityOutput = wModality' * dOmega{struct};
        end

        function gradientOutput = projectGradient(this,dij,struct,fGrad,w)
            dOmegaQt =  this.getSubQuantity([this.modality, 'dOmegaSqrtBeta']);
            dOmega = dOmegaQt.getResult(dij,w);
            
            gradientOutput = 2 * fGrad{struct} * dOmega{struct};
        end
    end
end