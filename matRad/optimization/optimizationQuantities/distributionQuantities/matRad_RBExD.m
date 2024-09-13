classdef matRad_RBExD < matRad_DistributionQuantity

    properties (Constant)
        quantityName = 'RBExD';
        requiredSubquantities = {'effect'};

    end

    methods
        function this = matRad_RBExD(dij)
            if nargin>0
                supArg = {dij};
            else
                supArg = {};
            end
            
            this@matRad_DistributionQuantity(supArg{:});
        end

        function quantityOutput = computeQuantity(this, dij, scen,w)
            effect = this.subQuantities{1}.getResult(dij,w);

            quantityOutput = zeros(dij.doseGrid.numOfVoxels,1);
            [ctScen,~] = ind2sub(size(dij.physicalDose),scen);
            quantityOutput(dij.ixDose{ctScen}) = sqrt((effect{scen}(dij.ixDose{ctScen})./dij.bx{ctScen}(dij.ixDose{ctScen}))+(dij.gamma{ctScen}(dij.ixDose{ctScen}).^2)) - dij.gamma{ctScen}(dij.ixDose{ctScen});
            
        end

        function gradientOutput = projectGradient(this,dij,scen,fGrad,w)
            
            currRBEvalue = this.getResult(dij,w);

            [ctScen,~] = ind2sub(size(dij.physicalDose),scen);
 
            scaledEffect = currRBEvalue{scen} + dij.gamma{ctScen};

            fGradTemp = zeros(dij.doseGrid.numOfVoxels,1);
            fGradTemp(dij.ixDose{ctScen}) = fGrad{scen}(dij.ixDose{ctScen}) ./ (2*dij.bx{ctScen}(dij.ixDose{ctScen}).*scaledEffect(dij.ixDose{ctScen}));

            fGradRBE = cell(numel(this.useScenarios),1);
            fGradRBE{scen} = fGradTemp; 
            
            gradientOutput = this.subQuantities{1}.projectGradient(dij,scen,fGradRBE,w);
        end
    end
end