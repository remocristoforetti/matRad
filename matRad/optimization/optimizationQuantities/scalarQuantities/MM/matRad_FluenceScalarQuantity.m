classdef (Abstract) matRad_FluenceScalarQuantity < matRad_ScalarQuantity

    properties (Abstract, Constant)
        dijField;
    end

    methods

        function this = matRad_FluenceScalarQuantity(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});
        end

        function quantityOutput = computeQuantity(this, dij,struct,w)
            w   = w.(this.modality);
            quantityOutput = dij.(this.modality).(this.dijField{1}){struct} * w;
        end

        function gradientOutput = projectGradient(this,dij,struct,fGrad,~)
            % This will not work but should also be not needed
            gradientOutput = fGrad{struct} * dij.(this.modality).(this.dijField{1}){struct}';
        end
    end
end