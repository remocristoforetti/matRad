classdef (Abstract) matRad_OptimizationQuantity < handle

    properties (Abstract,Constant)
        quantityName;
        requiredSubquantities;
    end

    properties
        subQuantities;
        d;
        wGrad;
        wCache;
        wGradCache;
        scenarios;
    end

    methods
        function this = matRad_OptimizationQuantity()

        end

        function quantityOutput = getResult(this,dij,w)
            if ~isequal(this.wCache,w)
                this.d = arrayfun(@(scen) this.computeQuantity(dij,scen,w), this.scenarios, 'UniformOutput',false);
                this.wCache = w;
            end
            quantityOutput = this.d;
        end

        function gradOutput = getProjectedGradient(this,dij,fGrad,w)
            if ~isequal(this.wGradCache,w)
                this.wGrad = arrayfun(@(scen) this.projectGradient(dij,scen,fGrad,w), this.scenarios, 'UniformOutput',false);
                this.wGradCache = w;
            end
            gradOutput = this.wGrad;
        end

        function output = computeQuantity(~)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('Function needs to be implemented by subclass');
            output = [];
        end

         function output = projectGradient(~)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('Function needs to be implemented by subclass');
            output = [];
        end
    end

    methods (Static)
        function optiFunc = setBiologicalDosePrescriptions(optiFunc,alphaX,betaX)

        end
    end
end