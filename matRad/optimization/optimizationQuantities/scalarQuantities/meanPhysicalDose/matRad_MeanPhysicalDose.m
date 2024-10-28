classdef matRad_MeanPhysicalDose < matRad_ScalarQuantity

    properties (Constant)
        quantityName = 'meanPhysicalDose';
        requiredSubquantities = {};

    end

    methods
        function this = matRad_MeanPhysicalDose(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});


        end

        function quantityOutput = computeQuantity(this, dij, struct,w)
            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);

            quantityOutput = (1/N) * dij.physicalDoseJ{struct}' * w;
        end

        function gradientOutput = projectGradient(this,dij,struct,fGrad,~)

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);

            
            gradientOutput = (1/N)*fGrad{struct} *  dij.physicalDoseJ{struct};
        end

        function constJacobianOutput = projectConstraintJacobian(this,dij,struct,fJacob,~)
            
            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);

            constJacobianOutput = (1/N) * fJacob{struct} *  dij.physicalDoseJ{struct}';
        end
    end  
end