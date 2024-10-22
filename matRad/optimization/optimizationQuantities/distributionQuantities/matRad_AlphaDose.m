classdef matRad_AlphaDose < matRad_DistributionQuantity

    properties (Constant)
        quantityName = 'AlphaDose';
        requiredSubquantities = {};
    end

    properties

    end

    methods
        function this = matRad_AlphaDose(dij)
            if nargin>0
                supArg = {dij};
            else
                supArg = {};
            end
            
            this@matRad_DistributionQuantity(supArg{:});
        end

        function  quantityOutput = computeQuantity(~, dij, scen,w)
            quantityOutput = dij.mAlphaDose{scen}*w;
        end

        function  gradientProjectionOutput = projectGradient(~, dij, scen,fGrad,~)
                % This is used to project the gradient and use the fGrad to
                % sum over all the voxels for example.              
                gradientProjectionOutput = fGrad{scen}' * dij.mAlphaDose{scen};

        end

        function constJacobianOutput = projectConstraintJacobian(~,dij,fJacob,~)
            constJacobianOutput = fJacob{1}' * dij.mAlphaDose{1};
        end


    end
end