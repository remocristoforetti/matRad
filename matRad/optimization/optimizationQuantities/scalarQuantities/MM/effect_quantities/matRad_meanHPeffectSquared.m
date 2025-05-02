classdef (Abstract) matRad_meanHPeffectSquared < matRad_ScalarQuantity
    % This might be just a helper function for later on to compute the
    % squared terms appering in ther definitiopn of theh approximated
    % robust effect. MIgt not be necessary for now.
    properties
        
    end

    methods
        function this = matRad_meanHPeffectSquared(cst)
             if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end

            this@matRad_ScalarQuantity(supArg{:});
        end


        function quantityOutput = computeQuantity(this, dij,struct,w)
            
            effectSubQt = this.getSubQuantity(this.effectSubquantity);

            % This sshould output teh value for a "single scenario". There
            % are no scenarios here. This should just be a total value
            % computed over the scenarios
            effectDistribution = effectSubQt.getResult(dij,w);

            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);

            quantityOutput = (1/N^2)*effectDistribution{1}(currIdx)'*effectDistribution{1}(currIdx);

        end

        function gradientOutput = projectGradient(this,dij,struct,fGrad,w)

        end
    end
end