classdef matRad_ConstRBEDose < matRad_DijDistributionQuantity

    properties (Constant)
        quantityName = 'constRBEDose';
        requiredSubquantities = {};
        
        dijField = {'physicalDose'};
    end

    methods

        function this = matRad_ConstRBEDose(dij)

            if nargin>0
                supArg = {dij};
            else
                supArg = {};
            end
            
            this@matRad_DijDistributionQuantity(supArg{:});
        end

        function quantityOutput = computeQuantity(this, dij, scen,w)
            physicalDose = computeQuantity@matRad_DijDistributionQuantity(this,dij,scen,w);
            quantityOutput = physicalDose*dij.RBE;
        
        end

        function gradientOutput = projectGradient(this,dij,scen,fGrad,w)
            
            physicalDoseGradient = projectGradient@matRad_DijDistributionQuantity(this,dij,scen,fGrad,w);
            gradientOutput = physicalDoseGradient*dij.RBE;

        end


        function constJacobianOutput = projectConstraintJacobian(this,dij,fJacob,w)

            physicalDoseJacobian = projectConstraintJacobian@matRad_DijDistributionQuantity(this,dij,fJacob,w);
            constJacobianOutput = physicalDoseJacobian*dij.RBE;
        end


    end
end