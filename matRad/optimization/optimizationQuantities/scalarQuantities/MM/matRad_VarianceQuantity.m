classdef (Abstract) matRad_VarianceQuantity < matRad_ScalarQuantity

    properties (Abstract, Constant)
        vMeanQuantity
        meanQuantity
    end

    methods
        function this = matRad_VarianceQuantity(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});
        end

        function quantityOutput = computeQuantity(this, dij,struct,w)

            vMeanQt = this.getSubQuantity(this.vMeanQuantity);
            vMean   = vMeanQt.getResult(dij,w);

            meanQt    = this.getSubQuantity(this.meanQuantity);
            meanValue = meanQt.getResult(dij,w);

            quantityOutput = vMean{struct} - meanValue{struct}^2;
        end

        function gradientOutput = projectGradient(this,dij,struct,fGrad,w)

            vMeanQt = this.getSubQuantity(this.vMeanQuantity);
            %vMean   = vMeanQt.getResult(dij,w);

            meanQt    = this.getSubQuantity(this.meanQuantity);
            meanValue = meanQt.getResult(dij,w);

            % This is the gradient of the wOmegaw part
            vMeanGradient = vMeanQt.projectGradient(dij,struct,fGrad,w);

            % This is the fGrad for the meanValue part
            meanValuefGrad{struct} = (2 * fGrad{struct}) * meanValue{struct};
            
            % This is the gradient for the meanValue part
            meanValueGrad = meanQt.projectGradient(dij,struct,meanValuefGrad,w);

            gradientOutput = vMeanGradient - meanValueGrad;
        end

    end
end