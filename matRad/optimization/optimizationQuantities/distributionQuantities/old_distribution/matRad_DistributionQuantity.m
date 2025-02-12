classdef (Abstract) matRad_DistributionQuantity < matRad_OptimizationQuantity

    properties
        useScenarios;
    end

    methods
        function this = matRad_DistributionQuantity(dij)
            this@matRad_OptimizationQuantity();
            if exist('dij', 'var')
                this.initializeProperties(dij);
            end

        end

        function quantityOutput = getResult(this,dij,w)
            if ~isequal(this.wCache,w)
                this.d(this.useScenarios) = arrayfun(@(scen) this.computeQuantity(dij,scen,w), this.useScenarios, 'UniformOutput',false);
                this.wCache = w;
            end
            quantityOutput = this.d;
        end

        function gradOutput = getProjectedGradient(this,dij,fGrad,w)
            if ~isequal(this.wGradCache,w)
                this.wGrad(this.useScenarios) = arrayfun(@(scen) this.projectGradient(dij,scen,fGrad,w), this.useScenarios, 'UniformOutput',false);
                this.wGradCache = w;
            end
            gradOutput = this.wGrad;
        end

        % function constraintOutput = getConstraintResult(this,dij,w)
        %     if ~isequal(this.wConstraintCache,w)
        %         this.c = this.computeConstraint(dij,w);
        %         this.wConstraintCache = w;
        %     end
        %     constraintOutput = this.c;
        % end

        function constJacobianOutput = getProjectedJacobian(this,dij,fJacob,w)
            if ~isequal(this.wConstJacobianCache,w)
                this.wJacob{1} = this.projectConstraintJacobian(dij,fJacob,w);
                this.wConstJacobianCache = w;
            end
            constJacobianOutput = this.wJacob;
        end

        function initializeProperties(this,dij)
            
            % This is quite a mess
            if isfield(dij, 'radiationModalities')

                modalityName = dij.radiationModalities{1};
                if isfield(dij.(modalityName), 'physicalDose') && ~isempty(dij.(modalityName).physicalDose{1})
                    distributionQuantity = 'physicalDose';
                elseif isfield(dij.(modalityName), 'physicalDoseExp') && ~isempty(dij.(modalityName).physicalDoseExp{1})
                    distributionQuantity = 'physicalDoseExp';
                elseif isfield(dij.(modalityName), 'mAlphaDose') && ~isempty(dij.(modalityName).mAlphaDose{1})
                    distributionQuantity = 'mAlphaDose';
                elseif isfield(dij.(modalityName), 'mAlphaDoseExp') && ~isempty(dij.(modalityName).mAlphaDoseExp{1})
                    distributionQuantity = 'mAlphaDoseExp';
                end
 
                this.d     = cell(size(dij.(modalityName).(distributionQuantity)));
                this.wGrad = cell(size(dij.(modalityName).(distributionQuantity)));
                this.wJacob = cell(1);


            else
                if isfield(dij, 'physicalDose') && ~isempty(dij.physicalDose{1})
                    distributionQuantity = 'physicalDose';
                elseif isfield(dij, 'physicalDoseExp') && ~isempty(dij.physicalDoseExp{1})
                    distributionQuantity = 'physicalDoseExp';
                elseif isfield(dij, 'mAlphaDose') && ~isempty(dij.mAlphaDose{1})
                    distributionQuantity = 'mAlphaDose';
                elseif isfield(dij, 'mAlphaDoseExp') && ~isempty(dij.mAlphaDoseExp{1})
                    distributionQuantity = 'mAlphaDoseExp';
                end

                this.d     = cell(size(dij.(distributionQuantity)));
                this.wGrad = cell(size(dij.(distributionQuantity)));
                this.wJacob = cell(1);

            end

        end

    end
end