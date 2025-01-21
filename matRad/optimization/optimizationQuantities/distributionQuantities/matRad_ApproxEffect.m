classdef matRad_ApproxEffect < matRad_DistributionQuantity

    properties (Constant)
        quantityName = 'ApproxEffect';
        requiredSubquantities = {'AlphaDoseExp', 'SqrtBetaDoseExp'};
    end

    properties
        
    end

    methods
        function this = matRad_ApproxEffect(dij)
            if nargin>0
                supArg = {dij};
            else
                supArg = {};
            end
            
            this@matRad_DistributionQuantity(supArg{:});

        end

        function quantityOutput = computeQuantity(this, dij, scen,w)

            alphaDose = this.subQuantities{1}.getResult(dij,w);
            betaDose  = this.subQuantities{2}.getResult(dij,w);

            quantityOutput = alphaDose{scen} + betaDose{scen}.^2;
        end

        function gradientProjectionOutput = projectGradient(this,dij,scen,fGrad,w)
            

            AlphaDoseSubQuantity       = this.getSubQuantity('AlphaDoseExp');
            SqrtBetaDoseSubQuantity    = this.getSubQuantity('SqrtBetaDoseExp');


            alphaGrad = AlphaDoseSubQuantity.projectGradient(dij,scen,fGrad,w);
            betaDose  = SqrtBetaDoseSubQuantity.getResult(dij,w);
            
            fBetaGrad = cell(numel(this.useScenarios),1);
            fBetaGrad{scen} = 2 * (fGrad{scen} .* betaDose{scen});
            
            betaGrad = SqrtBetaDoseSubQuantity.projectGradient(dij,scen,fBetaGrad,w);

            vBias = alphaGrad;
            mPsi  = betaGrad;

            gradientProjectionOutput = (vBias + mPsi)';
        end

        function constJacobianOutput = projectConstraintJacobian(this,dij,fJacob,w)
            
            alphaSubQuantity    = this.getSubQuantity('AlphaDoseExp');
            sqrtBetaSubQuantity = this.getSubQuantity('SqrtBetaDoseExp');
            
            alphaJacob    = alphaSubQuantity.projectConstraintJacobian(dij,fJacob);
            sqrtBetaJacob = sqrtBetaSubQuantity.projectConstraintJacobian(dij,fJacob);
            sqrtBetaDose  = sqrtBetaSubQuantity.getResult(dij,w);

            fJacobBeta = {2 * sqrtBetaDose{1} .* fJacob{1}};
            betaJacob = sqrtBetaSubQuantity.projectConstraintJacobian(dij,fJacobBeta);
            
            constJacobianOutput = alphaJacob + betaJacob;
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