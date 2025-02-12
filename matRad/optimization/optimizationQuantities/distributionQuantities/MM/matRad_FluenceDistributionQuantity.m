classdef (Abstract) matRad_FluenceDistributionQuantity < matRad_DistributionQuantity

    properties (Abstract, Constant)
        dijField;
    end

    methods

        function this = matRad_FluenceDistributionQuantity(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_DistributionQuantity(supArg{:});
        end

        function quantityOutput = computeQuantity(this, dij,scen,w)
            w   = w.(this.modality);
            quantityOutput = dij.(this.modality).(this.dijField{1}){scen} * w;
        end

        function gradientOutput = projectGradient(this,dij,scen,fGrad,~)
            % This will not work but should also be not needed
            gradientOutput = (fGrad{scen}' * dij.(this.modality).(this.dijField{1}){scen})';
        end
    end
end