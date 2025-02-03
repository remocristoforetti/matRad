classdef matRad_SMtotalVarianceAlpha < matRad_ScalarQuantity

    properties (Constant)

    end

    methods
        function this = matRad_SMtotalVarianceAlpha(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});
        end

        function quantityOutput = computeQuantity(this, dij, struct,w)
            
            wModality = w.(this.modality);

            dOmegaAlphaQt =  this.getSubQuantity([this.modality, 'dOmegaAlpha']);
            dOmegaAlpha = dOmegaAlphaQt.getResult(dij,w);
            
            quantityOutput = wModality' * dOmegaAlpha{struct};
        end

        function gradientOutput = projectGradient(this,dij,struct,fGrad,w)
            dOmegaAlphaQt =  this.getSubQuantity([this.modality, 'dOmegaAlpha']);
            dOmegaAlpha = dOmegaAlphaQt.getResult(dij,w);
            
            gradientOutput = 2 * fGrad{struct} * dOmegaAlpha{struct};
        end
    end
end