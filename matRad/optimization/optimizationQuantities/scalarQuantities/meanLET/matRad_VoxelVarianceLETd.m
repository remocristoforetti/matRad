classdef matRad_VoxelVarianceLETd < matRad_ScalarQuantity

    properties (Constant)
        quantityName = 'VoxelVarianceLETd';
        requiredSubquantities = {'vLETd', 'MeanLETd'};

    end

    methods
        function this = matRad_VoxelVarianceLETd(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});


        end

        function quantityOutput = computeQuantity(this, dij, struct,w)
            
            MeanLETdSubQuantity     = this.getSubQuantity('MeanLETd');
            LETdOmegaSubQuantity    = this.getSubQuantity('vLETd');

            meanLETd  = MeanLETdSubQuantity.getResult(dij,w);
            LEtdOmega = LETdOmegaSubQuantity.getResult(dij,w);

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);

            quantityOutput = (1/N)*(LEtdOmega{struct}) - meanLETd{struct}^2;
        end

        function gradientOutput = projectGradient(this,dij,struct,fGrad,w)

            MeanLETdSubQuantity     = this.getSubQuantity('MeanLETd');
            LETdOmegaSubQuantity    = this.getSubQuantity('vLETd');

            meanLETd  = MeanLETdSubQuantity.getResult(dij,w);

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);

            LETdOmegaGrad = LETdOmegaSubQuantity.projectGradient(dij,struct,fGrad,w);

            meanLETdfGrad{struct} = (2 * fGrad{struct}) * meanLETd{struct};
            
            meanLETdGrad = MeanLETdSubQuantity.projectGradient(dij,struct,meanLETdfGrad,w);

            gradientOutput = (1/N)*(LETdOmegaGrad) - meanLETdGrad;
        end
    end

  
end