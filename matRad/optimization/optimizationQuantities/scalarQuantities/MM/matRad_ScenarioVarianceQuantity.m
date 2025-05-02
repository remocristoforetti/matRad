classdef (Abstract) matRad_ScenarioVarianceQuantity < matRad_ScalarQuantity

    properties (Abstract, Constant)
        vMeanQuantity
        meanDistributionQuantity
    end

    methods
        function this = matRad_ScenarioVarianceQuantity(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});
        end

        function quantityOutput = computeQuantity(this, dij,struct,w)

            % This is the (1/N) w * Omega *w part
            vMeanQt = this.getSubQuantity(this.vMeanQuantity);
            vMean   = vMeanQt.getResult(dij,w);

            % This is just the distribution quantity
            meanQt    = this.getSubQuantity(this.meanDistributionQuantity);
            meanValue = meanQt.getResult(dij,w);

            % Get voxels in the current structure
            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);

            quantityOutput = vMean{struct} - (1/N).*meanValue{1}(currIdx)' * meanValue{1}(currIdx);
        end

        function gradientOutput = projectGradient(this,dij,struct,fGrad,w)

            vMeanQt = this.getSubQuantity(this.vMeanQuantity);

            meanQt    = this.getSubQuantity(this.meanDistributionQuantity);
            meanValue = meanQt.getResult(dij,w);

            % Get voxels in the current structure
            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);

            % This is the gradient of the wOmegaw part
            vMeanGradient = vMeanQt.projectGradient(dij,struct,fGrad,w);

            % This is the fGrad for the meanValue part
            tmpStructMeanValue = zeros(size(meanValue{1}));
            tmpStructMeanValue(currIdx) = meanValue{1}(currIdx);

            meanValuefGrad{struct} = (2 *(1/N)*fGrad{struct}) * tmpStructMeanValue;
            
            % This is the gradient for the meanValue part
            meanValueGrad = meanQt.projectGradient(dij,struct,meanValuefGrad,w);

            gradientOutput = vMeanGradient - meanValueGrad;
        end

    end
end