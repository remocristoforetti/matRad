classdef matRad_Effect < matRad_OptimizationQuantity

    properties (Constant)
        quantityName = 'effect';
        requiredSubquantities = {'AlphaDose', 'SqrtBetaDose'};
    end

    properties
        
    end

    methods
        function this = matRad_Effect()
            this@matRad_OptimizationQuantity();
        end

        function quantityOutput = computeQuantity(this, dij, scen,w)

            alphaDose = this.subQuantities{1}.getResult(dij,w);
            betaDose  = this.subQuantities{2}.getResult(dij,w);

            quantityOutput = alphaDose{scen} + betaDose{scen}.^2;
        end

        function gradientOutput = projectGradient(this,dij,scen,fGrad,w)
            
            alphaGrad = this.subQuantities{1}.getProjectedGradient(dij,[],w);
            betaGrad  = this.subQuantities{2}.getProjectedGradient(dij,[],w);
            betaDose  = this.subQuantities{2}.getResult(dij,w);
            
            vBias = fGrad{scen}' * alphaGrad{scen};
            mPsi  = 2 *(fGrad{scen} .* betaDose{scen})' * betaGrad{scen};  
            
            gradientOutput = (vBias + mPsi)';
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