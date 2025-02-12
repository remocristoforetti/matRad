classdef (Abstract) matRad_vMeanScalarQuantity  < matRad_vTotScalarQuantity
    % This implements the average variance part (normalizes by N the vTotScalarQuantity)
    properties

    end

    methods
        function this = matRad_vMeanScalarQuantity(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_vTotScalarQuantity(supArg{:});            
        end

        function quantityOutput = computeQuantity(this, dij,struct,w)

             % Call superclass computation, this just performs
             % dij.(modality).(fieldName){struct} * w
             quantityOutput = computeQuantity@matRad_vTotScalarQuantity(this,dij,struct,w);

             % Get the number of voxels in this struct
             currIdx = cat(1,this.cst{struct,4}{:});
             currIdx = unique(currIdx);
             N = numel(currIdx);

             % Normalize
             quantityOutput = (1/N) * quantityOutput;
         end

         function gradientOutput = projectGradient(this,dij,struct,fGrad,w)
            
            gradientOutput = projectGradient@matRad_vTotScalarQuantity(this,dij,struct,fGrad,w);
            
            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);

            gradientOutput = (1/N) * gradientOutput;

        end
    end
end