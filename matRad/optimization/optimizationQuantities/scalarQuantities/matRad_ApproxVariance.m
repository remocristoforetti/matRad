classdef matRad_ApproxVariance < matRad_ScalarQuantity

    properties (Constant)
        quantityName = 'vTotApprox';
        requiredSubquantities = {'vAlpha', 'ApproxEffect'};

    end

    methods
        function this = matRad_ApproxVariance(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});


        end

        function quantityOutput = computeQuantity(this, dij, struct,w)
            
            vAlphaSubQuantity       = this.getSubQuantity('vAlpha');
            ApproxEffectSubquantity = this.getSubQuantity('ApproxEffect');
            
            vAlpha = vAlphaSubQuantity.getResult(dij,w);
            approxEffect = ApproxEffectSubquantity.getResult(dij,w);

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            
            quantityOutput = vAlpha{struct} - approxEffect{1}(currIdx)'*approxEffect{1}(currIdx);
        end

        function gradientOutput = projectGradient(this,dij,struct,fGrad,w)

            vAlphaSubQuantity       = this.getSubQuantity('vAlpha');
            ApproxEffectSubquantity = this.getSubQuantity('ApproxEffect');

            dOmegaAlpha = vAlphaSubQuantity.projectGradient(dij,struct,fGrad,w);

            approxEffect = ApproxEffectSubquantity.getResult(dij,w);

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            
            tmpApproxEffect = zeros(size(approxEffect{1}));
            tmpApproxEffect(currIdx) = approxEffect{1}(currIdx);
            approxEffect{1} = tmpApproxEffect;

            fGradEffect{1} = 2 * (fGrad{struct} .* approxEffect{1});

            dEffect     = ApproxEffectSubquantity.projectGradient(dij,1,fGradEffect,w);
           
            gradientOutput =dOmegaAlpha - dEffect; %cellfun(@(structdOmega) fGrad * structdOmega, dOmega, 'UniformOutput', false);
        end
    end
end