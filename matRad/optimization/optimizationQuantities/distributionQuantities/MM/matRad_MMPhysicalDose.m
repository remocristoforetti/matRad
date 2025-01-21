classdef matRad_MMPhysicalDose < matRad_DistributionQuantity

    properties
        SF;
    end

    properties (Constant)
        quantityName = 'MMPhysicalDose';
        requiredSubquantities = {'physicalDoseProtons', 'physicalDosePhotons'};
    end

    methods
        function this = matRad_MMPhysicalDose()
            this@matRad_DistributionQuantity();
        end

         function quantityOutput = computeQuantity(this, dij, scen,w)
            
            photonSubQuantity = this.getSubQuantity('physicalDosePhotons');
            protonSubQuantity = this.getSubQuantity('physicalDoseProtons');

            photonDose = photonSubQuantity.getResult(dij,w);
            protonDose = protonSubQuantity.getResult(dij,w);

            quantityOutput = photonDose{scen}*[this.SF.('photons')]' + protonDose{scen}*[this.SF.('protons')]';

         end

         function gradientOutput = projectGradient(this,dij,scen,fGrad,~)

            photonSubQuantity = this.getSubQuantity('physicalDosePhotons');
            protonSubQuantity = this.getSubQuantity('physicalDoseProtons');

            photonGradient = photonSubQuantity.projectGradient(dij,scen,fGrad);
            protonGradient = protonSubQuantity.projectGradient(dij,scen,fGrad);

            photonGradient = photonGradient.*[this.SF.photons]';
            protonGradient = protonGradient.*[this.SF.photons]';
            
            % modalities = dij.radiationModalities;
            % for modalityIdx=1:dij.numOfModalities
            %     modalityName = dij.radiationModalities

            % Pay attention to the order of these modalities! Should
            % correspond with the weight splitting order
            
            gradientOutput = [protonGradient(:); photonGradient(:)];
         end
    end
end