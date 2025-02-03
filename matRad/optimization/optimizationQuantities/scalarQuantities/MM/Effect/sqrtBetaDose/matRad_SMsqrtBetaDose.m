classdef matRad_SMsqrtBetaDose < matRad_ScalarQuantity
    
    properties (Constant)
        
    end

    methods
        function this = matRad_SMsqrtBetaDose(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});
        end
 
        function quantityOutput = computeQuantity(this, dij, struct,w)

            w   = w.(this.modality);

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);

            quantityOutput = (1/N) * dij.(this.modality).sqrtBetaDoseJ{struct}' * w;

        end

        function gradientOutput = projectGradient(this,dij,struct,fGrad,~)
            
            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);


            gradientOutput = (1/N)*fGrad{struct} *  dij.(this.modality).sqrtBetaDoseJ{struct};

        end
    end

end