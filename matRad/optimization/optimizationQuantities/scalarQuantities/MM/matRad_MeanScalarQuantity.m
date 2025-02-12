classdef (Abstract) matRad_MeanScalarQuantity < matRad_FluenceScalarQuantity

    properties

    end

    methods
        function this = matRad_MeanScalarQuantity(cst)
            if nargin>0
                supArg = {cst};
            else
                supArg = {};
            end
            
            this@matRad_FluenceScalarQuantity(supArg{:});            
        end

         function quantityOutput = computeQuantity(this, dij,struct,w)

             % Call superclass computation, this just performs
             % dij.(modality).(fieldName){struct} * w
             quantityOutput = computeQuantity@matRad_FluenceScalarQuantity(this,dij,struct,w);

             % Get the number of voxels in this struct
             currIdx = cat(1,this.cst{struct,4}{:});
             currIdx = unique(currIdx);
             N = numel(currIdx);

             % Normalize
             quantityOutput = (1/N) * quantityOutput;
         end

         function gradientOutput = projectGradient(this,dij,struct,fGrad,~)
            
            gradientOutput = projectGradient@matRad_FluenceScalarQuantity(this,dij,struct,fGrad);
            
            currIdx = cat(1,this.cst{struct,4}{:});
            currIdx = unique(currIdx);
            N = numel(currIdx);

            gradientOutput = (1/N) * gradientOutput;

        end

    end
end