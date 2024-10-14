classdef matRad_MeanLETd < matRad_ScalarQuantity

    properties (Constant)
        quantityName = 'MeanLETd';
        requiredSubquantities = {};

    end

    methods
        function this = matRad_MeanLETd(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_ScalarQuantity(supArg{:});


        end

        function quantityOutput = computeQuantity(this, dij, struct,w)
            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);

            quantityOutput = (1/N) * dij.mLETdJ{struct}' * w;
        end

        function gradientOutput = projectGradient(this,dij,struct,fGrad,~)

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);

            
            gradientOutput = (1/N)*fGrad{struct} *  dij.mLETdJ{struct};
        end
    end  
end