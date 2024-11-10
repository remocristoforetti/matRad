classdef matRad_SqrtBetaOnlyVariance < matRad_ScalarQuantity

    properties (Constant)
        quantityName = 'SqrtBetaOnlyVariance';
        requiredSubquantities = {'SqrtBetaDoseExp', 'vBeta'};

    end

    methods
        function this = matRad_SqrtBetaOnlyVariance(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});


        end

        function quantityOutput = computeQuantity(this, dij, struct,w)
            
            AlphaDoseSubQuantity     = this.getSubQuantity('SqrtBetaDoseExp');
            AlphaOmegaSubQuantity    = this.getSubQuantity('vBeta');

            alphaDose  = AlphaDoseSubQuantity.getResult(dij,w);
            alphaOmega = AlphaOmegaSubQuantity.getResult(dij,w);

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);

            quantityOutput = (1/N)*(alphaOmega{struct} - (1/N)*(sum(alphaDose{1}(currIdx))*sum(alphaDose{1}(currIdx))));
        end

        function gradientOutput = projectGradient(this,dij,struct,fGrad,w)

            AlphaDoseSubQuantity     = this.getSubQuantity('SqrtBetaDoseExp');
            AlphaOmegaSubQuantity    = this.getSubQuantity('vBeta');

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);

            alphaOmegaGrad = AlphaOmegaSubQuantity.projectGradient(dij,struct,fGrad,w);

            alphaDose = AlphaDoseSubQuantity.getResult(dij,w);

            fAlphaGrad{1} = zeros(dij.doseGrid.numOfVoxels,1);
            fAlphaGrad{1}(currIdx) = (2 * fGrad{struct})*(sum(alphaDose{1}(currIdx))/N);
            
            alphaDoseGrad = AlphaDoseSubQuantity.projectGradient(dij,1,fAlphaGrad,w)';

            gradientOutput = (1/N)*(alphaOmegaGrad - alphaDoseGrad);
        end
    end

  
end