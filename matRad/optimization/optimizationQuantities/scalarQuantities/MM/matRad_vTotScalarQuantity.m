classdef (Abstract) matRad_vTotScalarQuantity < matRad_ScalarQuantity

    properties

    end

    methods
        function this = matRad_vTotScalarQuantity(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});            
        end


        function quantityOutput = computeQuantity(this, dij,struct,w)

             % Get the Omega * w part
             dOmegaQt = this.getSubQuantity(this.dOmegaSubQt);
             dOmega = dOmegaQt.getResult(dij,w);
 
             % Select modality weights
             w   = w.(this.modality);
             
             % Compute quantity
             quantityOutput = w' * dOmega{struct};
         end

         function gradientOutput = projectGradient(this,dij,struct,fGrad,w)
            
             % Get the Omega * w part
            dOmegaQt = this.getSubQuantity(this.dOmegaSubQt);
            dOmega = dOmegaQt.getResult(dij,w);

            gradientOutput = 2 * fGrad{struct} * dOmega{struct};

         end

         function constJacobianOutput = projectConstraintJacobian(this,dij,struct,fJacob,w)
            
             % Get the Omega * w part
            dOmegaQt = this.getSubQuantity(this.dOmegaSubQt);
            dOmega = dOmegaQt.getResult(dij,w);

            constJacobianOutput = 2 * (dOmega{struct} * fJacob{struct})';

         end
    end
end