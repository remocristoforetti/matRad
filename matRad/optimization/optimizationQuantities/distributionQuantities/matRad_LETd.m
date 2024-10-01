classdef matRad_LETd < matRad_DistributionQuantity

    properties (Constant)
        quantityName = 'LETd';
        requiredSubquantities = {};
    end

    properties

    end

    methods
        function this = matRad_LETd(dij)
            if nargin>0
                supArg = {dij};
            else
                supArg = {};
            end
            
            this@matRad_DistributionQuantity(supArg{:});
        end

        function  quantityOutput = computeQuantity(~, dij, scen,w)
            quantityOutput = dij.mLETDose{scen}*w;
        end

        function  gradientProjectionOutput = projectGradient(~, dij, scen,fGrad,~)
                % This is used to project the gradient and use the fGrad to
                % sum over all the voxels for example.              
                gradientProjectionOutput = (fGrad{scen}' * dij.mLETDose{scen})';
        end
    end
end