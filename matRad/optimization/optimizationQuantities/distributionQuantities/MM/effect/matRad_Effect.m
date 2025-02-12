classdef (Abstract) matRad_Effect < matRad_DistributionQuantity

    properties (Abstract, Constant)
        alphaSubQuantity
        betaSubQuantity
    end

    methods
        function this = matRad_Effect(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_DistributionQuantity(supArg{:});
        end

         function quantityOutput = computeQuantity(this, dij, scen,w)
            
            alphaQuantity        = this.getSubQuantity(this.alphaSubQuantity);
            betaQuantity         = this.getSubQuantity(this.betaSubQuantity);

            alphaDose    = alphaQuantity.getResult(dij,w);
            betaDose     = betaQuantity.getResult(dij,w);

            quantityOutput = alphaDose{scen} + betaDose{scen}.^2;
         end

         function gradientOutput = projectGradient(this,dij,scen,fGrad,w)

            alphaQuantity        = this.getSubQuantity(this.alphaSubQuantity);
            betaQuantity         = this.getSubQuantity(this.betaSubQuantity);
            
            alphaGrad = alphaQuantity.projectGradient(dij,scen,fGrad,w);
            
            betaDose  = betaQuantity.getResult(dij,w);
            
            fBetaGrad = cell(numel(this.useScenarios),1);
            fBetaGrad{scen} = 2 * (fGrad{scen} .* betaDose{scen});

            betaGrad  = betaQuantity.projectGradient(dij,scen,fBetaGrad,w);

            gradientOutput = alphaGrad + betaGrad;
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