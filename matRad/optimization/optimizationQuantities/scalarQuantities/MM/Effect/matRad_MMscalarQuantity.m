classdef matRad_MMscalarQuantity < matRad_ScalarQuantity

    properties
        SF;
    end

    properties (Constant)
        %quantityName = 'MMmeanEffect';
        %requiredSubquantities = {'protonsMeanEffect', 'photonsMeanEffect'};
    
    end

    methods
        function this = matRad_MMscalarQuantity(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});
            
        end

        function quantityOutput = computeQuantity(this, dij, struct,w)

            protonQt = this.getSubQuantity(this.protonSubquantity);
            photonQt = this.getSubQuantity(this.photonSubquantity);

            protonResult = protonQt.getResult(dij,w);
            photonResult = photonQt.getResult(dij,w);

            quantityOutput = protonResult{struct}*[this.SF.('protons')]' + photonResult{struct}*[this.SF.('photons')]';
            
        end

        function gradientOutput = projectGradient(this,dij,struct,fGrad,w)

            protonQt = this.getSubQuantity(this.protonSubquantity);
            photonQt = this.getSubQuantity(this.photonSubquantity);
    
            protonGrad = protonQt.projectGradient(dij,struct,fGrad,w);
            photonGrad = photonQt.projectGradient(dij,struct,fGrad,w);

            protonGrad = protonGrad.*[this.SF.protons];
            photonGrad = photonGrad.*[this.SF.photons];
            
            % This is not robust agains modality order! Need to make this
            % consistent
            gradientOutput = [protonGrad(:); photonGrad(:)];

        end
    end

    methods
        
        function this = set.SF(this, value)

            this.SF = this.setSF(value);

        end
    end

    methods (Access = protected)

        function SF = setSF(this,value)
            % Some MM subquantities migh override this, to get squared
            % number of fractions for example
            SF = value;
        end
    end
end