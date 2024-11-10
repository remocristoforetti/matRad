classdef matRad_VoxelVarianceAlpha < matRad_ScalarQuantity

    properties (Constant)
        quantityName = 'VoxelVarianceAlpha';
        requiredSubquantities = {'vAlpha', 'MeanAlpha'};

    end

    methods
        function this = matRad_VoxelVarianceAlpha(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});


        end

        function quantityOutput = computeQuantity(this, dij, struct,w)
            
            MeanAlphaSubQuantity     = this.getSubQuantity('MeanAlpha');
            AlphaOmegaSubQuantity    = this.getSubQuantity('vAlpha');

            meanAlpha  = MeanAlphaSubQuantity.getResult(dij,w);
            alphaOmega = AlphaOmegaSubQuantity.getResult(dij,w);

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);

            quantityOutput = (1/N)*(alphaOmega{struct}) - meanAlpha{struct}^2;
        end

        function gradientOutput = projectGradient(this,dij,struct,fGrad,w)

            MeanAlphaSubQuantity     = this.getSubQuantity('MeanAlpha');
            AlphaOmegaSubQuantity    = this.getSubQuantity('vAlpha');

            meanAlpha  = MeanAlphaSubQuantity.getResult(dij,w);

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);

            alphaOmegaGrad = AlphaOmegaSubQuantity.projectGradient(dij,struct,fGrad,w);

            meanAlphafGrad{struct} = (2 * fGrad{struct}) * meanAlpha{struct};
            
            meanAlphaGrad = MeanAlphaSubQuantity.projectGradient(dij,struct,meanAlphafGrad,w);

            gradientOutput = (1/N)*(alphaOmegaGrad) - meanAlphaGrad;
        end
    end

  
end