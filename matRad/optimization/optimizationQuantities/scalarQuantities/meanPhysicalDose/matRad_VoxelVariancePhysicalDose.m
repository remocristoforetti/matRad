classdef matRad_VoxelVariancePhysicalDose < matRad_ScalarQuantity

    properties (Constant)
        quantityName = 'VoxelVariancePhysicalDose';
        requiredSubquantities = {'vTotReduced', 'meanPhysicalDose'};

    end

    methods
        function this = matRad_VoxelVariancePhysicalDose(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});


        end

        function quantityOutput = computeQuantity(this, dij, struct,w)
            
            meanPhysicalDoseSubQuantity     = this.getSubQuantity('meanPhysicalDose');
            physicalDoseVarianceSubQuantity    = this.getSubQuantity('vTotReduced');

            meanPhysicalDose     = meanPhysicalDoseSubQuantity.getResult(dij,w);
            physicalDoseVariance = physicalDoseVarianceSubQuantity.getResult(dij,w);

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);

            quantityOutput = (1/N)*(physicalDoseVariance{struct}) - meanPhysicalDose{struct}^2;
        end

        function gradientOutput = projectGradient(this,dij,struct,fGrad,w)

            meanPhysicalDoseSubQuantity     = this.getSubQuantity('meanPhysicalDose');
            physicalDoseVarianceSubQuantity = this.getSubQuantity('vTotReduced');

            meanPhysicalDose  = meanPhysicalDoseSubQuantity.getResult(dij,w);

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);

            physicalDoseVarianceGrad = physicalDoseVarianceSubQuantity.projectGradient(dij,struct,fGrad,w);

            meanPhysicalDosefGrad{struct} = (2 * fGrad{struct}) * meanPhysicalDose{struct};
            
            meanPhysicalDoseGrad = meanPhysicalDoseSubQuantity.projectGradient(dij,struct,meanPhysicalDosefGrad,w);

            gradientOutput = (1/N)*(physicalDoseVarianceGrad) - meanPhysicalDoseGrad;
        end

        function constJacobianOutput = projectConstraintJacobian(this,dij,struct,fJacob,w)
            
            meanPhysicalDoseSubQuantity     = this.getSubQuantity('meanPhysicalDose');
            physicalDoseVarianceSubQuantity = this.getSubQuantity('vTotReduced');

            meanPhysicalDose  = meanPhysicalDoseSubQuantity.getResult(dij,w);

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);

            physicalDoseVarianceGrad = physicalDoseVarianceSubQuantity.projectConstraintJacobian(dij,struct,fJacob,w);

            meanPhysicalDosefGrad{struct} = (2 * fJacob{struct}) * meanPhysicalDose{struct};
            
            meanPhysicalDoseGrad = meanPhysicalDoseSubQuantity.projectConstraintJacobian(dij,struct,meanPhysicalDosefGrad,w);

            constJacobianOutput = (1/N)*(physicalDoseVarianceGrad) - meanPhysicalDoseGrad;
        end
    end

  
end