classdef (Abstract) matRad_ScalarQuantity < matRad_OptimizationQuantity

    properties
        useStructs;
    end

    methods
        function this = matRad_ScalarQuantity(cst)
            this@matRad_OptimizationQuantity();
            if exist('cst', 'var')
                this.initializeProperties(cst);
            end
        end

        function quantityOutput = getResult(this,dij,w)
            if ~isequal(this.wCache,w)
                this.d(this.useStructs) = arrayfun(@(struct) this.computeQuantity(dij,struct,w), this.useStructs, 'UniformOutput',false);
                this.wCache = w;
            end
            quantityOutput = this.d;
        end

        function gradOutput = getProjectedGradient(this,dij,fGrad,w)
            if ~isequal(this.wGradCache,w)
                this.wGrad(this.useStructs) = arrayfun(@(struct) this.projectGradient(dij,struct,fGrad,w), this.useStructs, 'UniformOutput',false);
                this.wGradCache = w;
            end
            gradOutput = this.wGrad;
        end

        function initializeProperties(this,cst)
            this.d     = cell(size(cst),1);
            this.wGrad = cell(size(cst),1);
        end

    end
end