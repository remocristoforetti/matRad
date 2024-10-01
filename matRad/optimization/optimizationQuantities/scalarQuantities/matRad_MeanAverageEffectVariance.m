classdef matRad_MeanAverageEffectVariance < matRad_ScalarQuantity

    properties (Constant)
        quantityName = 'MeanAverageEffectVariance';
        requiredSubquantities = {'MeanAverageEffect', 'vAlpha'};

    end

    methods
        function this = matRad_MeanAverageEffectVariance(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});


        end

        function quantityOutput = computeQuantity(this, dij, struct,w)
            
            MAESubQuantity           = this.getSubQuantity('MeanAverageEffect');
            AlphaOmegaSubQuantity    = this.getSubQuantity('vAlpha');

            meanAverageEffect = MAESubQuantity.getResult(dij,w);
            alphaOmega        = AlphaOmegaSubQuantity.getResult(dij,w);

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);

            quantityOutput = (1/numel(currIdx))*(alphaOmega{struct}) - meanAverageEffect{struct}^2;
        end

        function gradientOutput = projectGradient(this,dij,struct,fGrad,w)

            MAESubQuantity           = this.getSubQuantity('MeanAverageEffect');
            AlphaOmegaSubQuantity    = this.getSubQuantity('vAlpha');

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);

            meanAverageEffect = MAESubQuantity.getResult(dij,w);
            alphaGrad = AlphaOmegaSubQuantity.projectGradient(dij,struct,fGrad,w);

            fMAEGrad = cell(size(fGrad));
            fMAEGrad{struct} = 2 * meanAverageEffect{struct} * fGrad{struct};

            maeGrad = MAESubQuantity.projectGradient(dij,struct,fMAEGrad,w);

            gradientOutput = (1/numel(currIdx))*alphaGrad - maeGrad;
        end
    end

  
end