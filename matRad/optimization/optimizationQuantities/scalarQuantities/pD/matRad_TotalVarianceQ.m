classdef matRad_TotalVarianceQ < matRad_ScalarQuantity

    properties (Constant)
        quantityName = 'vTot';
        requiredSubquantities = {'dOmega'};

    end

    methods
        function this = matRad_TotalVarianceQ(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});


        end

        function quantityOutput = computeQuantity(this, dij, struct,w)
            
            dOmega = this.subQuantities{1}.getResult(dij,w);
            quantityOutput = w' * dOmega{struct};%cellfun(@(structdOmega) w' * structdOmega, dOmega, 'UniformOutput',false);
        end

        function gradientOutput = projectGradient(this,dij,struct,fGrad,w)
            dOmega = this.subQuantities{1}.getResult(dij,w);
            gradientOutput = 2 * fGrad{struct} * dOmega{struct}; %cellfun(@(structdOmega) fGrad * structdOmega, dOmega, 'UniformOutput', false);
        end

        function constJacobianOutput = projectConstraintJacobian(this,dij,struct,fJacob,w)
            dOmega = this.subQuantities{1}.getResult(dij,w);
            constJacobianOutput = 2 * (dOmega{struct} * fJacob{struct})';
        end

    end
end