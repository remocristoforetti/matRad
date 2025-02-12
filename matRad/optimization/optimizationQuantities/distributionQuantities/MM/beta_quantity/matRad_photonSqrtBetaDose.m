classdef matRad_photonSqrtBetaDose < matRad_DistributionQuantity

    properties (Constant)
        quantityName = 'photonsSqrtBetaDose';
        requiredSubquantities = {'photonsDose'};

        modality = 'photons';
    end

    methods
        function this = matRad_photonSqrtBetaDose(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_DistributionQuantity(supArg{:});
        end

        function quantityOutput = computeQuantity(this, dij,scen,w)
            
            pDQt = this.getSubQuantity('photonsDose');
            pD   = pDQt.getResult(dij,w);

            quantityOutput = sqrt(dij.(this.modality).bx{1}) .* pD{scen};
        end

        function gradientOutput = projectGradient(this,dij,scen,fGrad,~)
            
            pDQt = this.getSubQuantity('photonsDose');
            fBetaGrad{scen} = sqrt(dij.(this.modality).bx{1}) .* fGrad{1};
            pDGrad   = pDQt.projectGradient(dij,scen, fBetaGrad);

            gradientOutput = pDGrad;
        end
    end
end