classdef matRad_VoxelVarianceSqrtBeta < matRad_ScalarQuantity

    properties (Constant)
        quantityName = 'VoxelVarianceSqrtBeta';
        requiredSubquantities = {'vBeta', 'MeanSqrtBeta'};

    end

    methods
        function this = matRad_VoxelVarianceSqrtBeta(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});


        end

        function quantityOutput = computeQuantity(this, dij, struct,w)
            
            MeanSqrtBetaSubQuantity     = this.getSubQuantity('MeanSqrtBeta');
            BetaOmegaSubQuantity    = this.getSubQuantity('vBeta');

            meanSqrtBeta  = MeanSqrtBetaSubQuantity.getResult(dij,w);
            betaOmega = BetaOmegaSubQuantity.getResult(dij,w);

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);

            quantityOutput = (1/N)*(betaOmega{struct}) - meanSqrtBeta{struct}^2;
        end

        function gradientOutput = projectGradient(this,dij,struct,fGrad,w)

            MeanSqrtBetaSubQuantity     = this.getSubQuantity('MeanSqrtBeta');
            BetaOmegaSubQuantity    = this.getSubQuantity('vBeta');

            meanSqrtBeta  = MeanSqrtBetaSubQuantity.getResult(dij,w);

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);

            betaOmegaGrad = BetaOmegaSubQuantity.projectGradient(dij,struct,fGrad,w);

            meanSqrtBetafGrad{struct} = (2 * fGrad{struct}) * meanSqrtBeta{struct};
            
            meanSqrtBetaGrad = MeanSqrtBetaSubQuantity.projectGradient(dij,struct,meanSqrtBetafGrad,w);

            gradientOutput = (1/N)*(betaOmegaGrad) - meanSqrtBetaGrad;
        end
    
    end
end