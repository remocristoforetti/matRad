classdef matRad_SMvoxelVarianceAlphaDose < matRad_ScalarQuantity

    properties (Constant)

    end

    methods
        function this = matRad_SMvoxelVarianceAlphaDose(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});


        end

        function quantityOutput = computeQuantity(this, dij, struct,w)
                        
            meanAlphaDoseSubQuantity        = this.getSubQuantity([this.modality, 'AlphaDose']);
            alphaDoseVarianceSubQuantity    = this.getSubQuantity([this.modality, 'vTotAlpha']);

            meanAlphaDose     = meanAlphaDoseSubQuantity.getResult(dij,w);
            alphaDoseVariance = alphaDoseVarianceSubQuantity.getResult(dij,w);

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);

            quantityOutput = (1/N)*(alphaDoseVariance{struct}) - meanAlphaDose{struct}^2;
        end

        function gradientOutput = projectGradient(this,dij,struct,fGrad,w)

            meanAlphaDoseSubQuantity        = this.getSubQuantity([this.modality, 'AlphaDose']);
            alphaDoseVarianceSubQuantity    = this.getSubQuantity([this.modality, 'vTotAlpha']);

            meanAlphaDose     = meanAlphaDoseSubQuantity.getResult(dij,w);

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);

            alphaDoseVarianceGrad = alphaDoseVarianceSubQuantity.projectGradient(dij,struct,fGrad,w);

            meanAlphaDosefGrad{struct} = (2 * fGrad{struct}) * meanAlphaDose{struct};
            
            meanAlphaDoseGrad = meanAlphaDoseSubQuantity.projectGradient(dij,struct,meanAlphaDosefGrad,w);

            gradientOutput = (1/N)*(alphaDoseVarianceGrad) - meanAlphaDoseGrad;
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