classdef matRad_CovarianceAlphaSqrtBeta < matRad_ScalarQuantity

    properties (Constant)
        quantityName = 'CovAlphaSqrtBeta';
        requiredSubquantities = {'vAlphaSqrtBeta', 'MeanAlpha', 'MeanSqrtBeta'};

    end

    methods
        function this = matRad_CovarianceAlphaSqrtBeta(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});


        end

        function quantityOutput = computeQuantity(this, dij, struct,w)
            
            AlphaSqrtBetaOmegaSubQuantity = this.getSubQuantity('vAlphaSqrtBeta');
            MeanAlphaSubQuantity          = this.getSubQuantity('MeanAlpha');
            MeanSqrtBetaSubQuantity       = this.getSubQuantity('MeanSqrtBeta');

            omegaAlphaSqrtBeta = AlphaSqrtBetaOmegaSubQuantity.getResult(dij,w);
            meanAlphaDose      = MeanAlphaSubQuantity.getResult(dij,w);
            meanBetaDose       = MeanSqrtBetaSubQuantity.getResult(dij,w);

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);


            quantityOutput = (1/N) * (omegaAlphaSqrtBeta{struct}) - (meanAlphaDose{struct} * meanBetaDose{struct});
        end

        function gradientOutput = projectGradient(this,dij,struct,fGrad,w)

            AlphaSqrtBetaOmegaSubQuantity = this.getSubQuantity('vAlphaSqrtBeta');
            MeanAlphaSubQuantity          = this.getSubQuantity('MeanAlpha');
            MeanSqrtBetaSubQuantity       = this.getSubQuantity('MeanSqrtBeta');

            meanAlphaDose = MeanAlphaSubQuantity.getResult(dij,w);
            meanBetaDose  = MeanSqrtBetaSubQuantity.getResult(dij,w);

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);
            
            fGradTimeBeta{struct}  = fGrad{struct} * meanBetaDose{struct};
            fGradTimeAlpha{struct} = fGrad{struct} * meanAlphaDose{struct};
            

            meanAlphaGrad = MeanAlphaSubQuantity.projectGradient(dij,struct,fGradTimeBeta,w);
            meanBetaGrad  = MeanSqrtBetaSubQuantity.projectGradient(dij,struct,fGradTimeAlpha,w);
            
            alphaSqrtBetaOmegaGrad = AlphaSqrtBetaOmegaSubQuantity.projectGradient(dij,struct,fGrad,w);

            gradientOutput = (1/N)*(alphaSqrtBetaOmegaGrad) - (meanAlphaGrad + meanBetaGrad);
        end
    end

end