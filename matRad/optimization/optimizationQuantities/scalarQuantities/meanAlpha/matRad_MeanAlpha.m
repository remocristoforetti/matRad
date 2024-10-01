classdef matRad_MeanAlpha < matRad_ScalarQuantity

    properties (Constant)
        quantityName = 'MeanAlpha';
        requiredSubquantities = {};

    end

    methods
        function this = matRad_MeanAlpha(cst)
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

            quantityOutput = (1/N) * dij.alphaJExp{struct}' * w;
        end

        function gradientOutput = projectGradient(this,dij,struct,fGrad,~)

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);

            
            gradientOutput = (1/N)*fGrad{struct} *  dij.alphaJExp{struct};
        end
    end  
end