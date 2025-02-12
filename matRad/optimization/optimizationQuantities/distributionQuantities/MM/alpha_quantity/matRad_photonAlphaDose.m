classdef matRad_photonAlphaDose < matRad_DistributionQuantity

    properties (Constant)
        quantityName = 'photonsAlphaDose';
         requiredSubquantities = {'photonsDose'};

        modality = 'photons';
    end

    methods
        function this = matRad_photonAlphaDose(cst)
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

            quantityOutput = dij.(this.modality).ax{1} .* pD{scen};
        end

        function gradientOutput = projectGradient(this,dij,scen,fGrad,~)
            
            pDQt = this.getSubQuantity('photonsDose');
            fGradAlpha{scen} = dij.(this.modality).ax{1} .* fGrad{1};
            pDGrad   = pDQt.projectGradient(dij,scen,fGradAlpha);

            gradientOutput = pDGrad;
        end
    end
end