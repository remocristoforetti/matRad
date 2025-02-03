classdef matRad_PhysicalDoseSingleModality < matRad_DistributionQuantity

    properties (Constant)

    end

    methods
        function this = matRad_PhysicalDoseSingleModality(dij)

            if nargin>0
                supArg = {dij};
            else
                supArg = {};
            end
            
            this@matRad_DistributionQuantity(supArg{:});

        end

        function quantityOutput = computeQuantity(this, dij, scen,w)

            w   = w.(this.modality);

            if ~isempty(dij.(this.modality).physicalDose{scen})

                quantityOutput = dij.(this.modality).physicalDose{scen}*w;

            else

                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
                quantityOutput = [];
            end
        end

        function gradientOutput = projectGradient(this,dij,scen,fGrad,~)
            if ~isempty(dij.(this.modality).physicalDose{scen})
                gradientOutput = (fGrad{scen}' * dij.(this.modality).physicalDose{scen})';
            else
                gradientOutput = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end
        end

        function constJacobianOutput = projectConstraintJacobian(~,dij,fJacob,~)
            constJacobianOutput = fJacob{1}' * dij.physicalDose{1};
        end

    end

end