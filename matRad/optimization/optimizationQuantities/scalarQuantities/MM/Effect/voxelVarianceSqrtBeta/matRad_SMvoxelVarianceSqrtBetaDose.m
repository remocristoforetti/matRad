classdef matRad_SMvoxelVarianceSqrtBetaDose < matRad_ScalarQuantity

    properties (Constant)

    end

    methods
        function this = matRad_SMvoxelVarianceSqrtBetaDose(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});


        end

        function quantityOutput = computeQuantity(this, dij, struct,w)
                        
            meanSubQuantity        = this.getSubQuantity([this.modality, 'SqrtBetaDose']);
            varianceSubQuantity    = this.getSubQuantity([this.modality, 'vTotSqrtBeta']);

            meanDose     = meanSubQuantity.getResult(dij,w);
            doseVariance = varianceSubQuantity.getResult(dij,w);

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);

            quantityOutput = (1/N)*(doseVariance{struct}) - meanDose{struct}^2;
        end

        function gradientOutput = projectGradient(this,dij,struct,fGrad,w)

            meanSubQuantity        = this.getSubQuantity([this.modality, 'SqrtBetaDose']);
            varianceSubQuantity    = this.getSubQuantity([this.modality, 'vTotSqrtBeta']);

            meanDose     = meanSubQuantity.getResult(dij,w);

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);

            varianceGrad = varianceSubQuantity.projectGradient(dij,struct,fGrad,w);

            meanfGrad{struct} = (2 * fGrad{struct}) * meanDose{struct};
            
            meanGrad = meanSubQuantity.projectGradient(dij,struct,meanfGrad,w);

            gradientOutput = (1/N)*(varianceGrad) - meanGrad;
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