classdef matRad_MeanEffect < matRad_ScalarQuantity

    properties (Constant)
        quantityName = 'MeanEffect';
        requiredSubquantities = {'MeanAlpha', 'vBeta'};

    end

    methods
        function this = matRad_MeanEffect(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});


        end

        function quantityOutput = computeQuantity(this, dij, struct,w)
            
            MeanAlphaSubQuantity        = this.getSubQuantity('MeanAlpha');
            SqrtBetaOmegaSubQuantity    = this.getSubQuantity('vBeta');

            meanAlphaDose = MeanAlphaSubQuantity.getResult(dij,w);
            betaOmega = SqrtBetaOmegaSubQuantity.getResult(dij,w);

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);

            quantityOutput = meanAlphaDose{struct} + (1/N) * betaOmega{struct};
        end

        function gradientOutput = projectGradient(this,dij,struct,fGrad,w)

            MeanAlphaSubQuantity        = this.getSubQuantity('MeanAlpha');
            SqrtBetaOmegaSubQuantity    = this.getSubQuantity('vBeta');

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);
            

            meanAlphaGrad = MeanAlphaSubQuantity.projectGradient(dij,struct,fGrad,w);
            betaOmegaGrad = SqrtBetaOmegaSubQuantity.projectGradient(dij,struct,fGrad,w);

            gradientOutput = meanAlphaGrad + (1/N) * betaOmegaGrad;
        end
    end

    methods (Static)
        function optiFunc = setBiologicalDosePrescriptions(optiFunc,alphaX,betaX)
            doses = optiFunc.getDoseParameters();
            effect = alphaX*doses + betaX*doses.^2;
            optiFunc = optiFunc.setDoseParameters(effect);
        end
    end
end