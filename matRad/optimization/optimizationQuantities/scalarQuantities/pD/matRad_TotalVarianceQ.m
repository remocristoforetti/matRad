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
            gradientOutput = fGrad * dOmega{struct}; %cellfun(@(structdOmega) fGrad * structdOmega, dOmega, 'UniformOutput', false);
        end
    end
end