classdef matRad_MeanAverageEffect < matRad_ScalarQuantity

    properties (Constant)
        quantityName = 'MeanAverageEffect';
        requiredSubquantities = {'AlphaDoseExp', 'vBeta'};

    end

    methods
        function this = matRad_MeanAverageEffect(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});


        end

        function quantityOutput = computeQuantity(this, dij, struct,w)
            
            AlphaDoseSubQuantity        = this.getSubQuantity('AlphaDoseExp');
            SqrtBetaOmegaSubQuantity    = this.getSubQuantity('vBeta');

            alphaDose = AlphaDoseSubQuantity.getResult(dij,w);
            betaOmega = SqrtBetaOmegaSubQuantity.getResult(dij,w);

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);

            quantityOutput = (1/numel(currIdx))*(sum(alphaDose{1}(currIdx)) + betaOmega{struct});
        end

        function gradientOutput = projectGradient(this,dij,struct,fGrad,w)

            AlphaDoseSubQuantity        = this.getSubQuantity('AlphaDoseExp');
            SqrtBetaOmegaSubQuantity    = this.getSubQuantity('vBeta');

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);

            fAlphaGrad{1} = zeros(dij.doseGrid.numOfVoxels,1);
            fAlphaGrad{1}(currIdx) = fGrad{struct}/numel(currIdx);

            alphaDoseGrad = AlphaDoseSubQuantity.projectGradient(dij,1,fAlphaGrad,w)';
            betaOmegaGrad = SqrtBetaOmegaSubQuantity.projectGradient(dij,struct,fGrad,w);

            gradientOutput = alphaDoseGrad + (betaOmegaGrad)/numel(currIdx); %cellfun(@(structdOmega) fGrad * structdOmega, dOmega, 'UniformOutput', false);
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