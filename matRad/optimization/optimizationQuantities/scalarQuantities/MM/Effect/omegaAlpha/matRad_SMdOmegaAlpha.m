classdef matRad_SMdOmegaAlpha < matRad_ScalarQuantity

    properties (Constant)

    end

    methods

        function this = matRad_SMdOmegaAlpha()
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});
        end

        function quantityOutput = computeQuantity(this, dij,struct,w)
            w   = w.(this.modality);
            quantityOutput = dij.(this.modality).mAlphaDoseOmega{struct} * w; %cellfun(@(structOmega) structOmega*w, dij. 'UniformOutput',false);
        end

        function gradientOutput = projectGradient(this,dij,struct,fGrad,~)
            % This will not work but should also be not needed
            gradientOutput = dij.(this.modality).mAlphaDoseOmega{struct};
        end

    
    end
end