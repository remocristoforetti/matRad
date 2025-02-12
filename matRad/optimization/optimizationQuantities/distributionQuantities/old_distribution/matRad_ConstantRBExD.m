classdef matRad_ConstantRBExD < matRad_DistributionQuantity

    properties (Constant)
        quantityName = 'constantRBExD';
        requiredSubquantities = {};
    end

    methods
        function this = matRad_ConstantRBExD(dij)

            if nargin>0
                supArg = {dij};
            else
                supArg = {};
            end
            
            this@matRad_DistributionQuantity(supArg{:});

        end

        function quantityOutput = computeQuantity(~, dij, scen,w)
            if ~isempty(dij.physicalDose{scen})
                
                quantityOutput = dij.physicalDose{scen}*w*dij.RBE;
            
            else

                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
                quantityOutput = [];
            end
        end

        function gradientOutput = projectGradient(~,dij,scen,fGrad,~)
            if ~isempty(dij.physicalDose{scen})
                gradientOutput = (fGrad{scen}' * dij.physicalDose{scen}*dij.RBE)';
            else
                gradientOutput = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end
        end

        % function constraintOutput = computeConstraint()
        % 
        % end

        function constJacobianOutput = projectConstraintJacobian(~,dij,fJacob,~)
            constJacobianOutput = fJacob{1}' * dij.physicalDose{1}*dij.RBE;
        end
    end

end