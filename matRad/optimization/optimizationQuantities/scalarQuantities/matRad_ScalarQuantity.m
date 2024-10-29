classdef (Abstract) matRad_ScalarQuantity < matRad_OptimizationQuantity

    properties
        useStructsOptimization;
        useStructsConstraint;
        cst;
    end

    methods
        function this = matRad_ScalarQuantity(cst)
            this@matRad_OptimizationQuantity();
            if exist('cst', 'var')
                this.initializeProperties(cst);
            end
        end

        function quantityOutput = getResult(this,dij,w)
            if ~isequal(this.wCache,w)
                % For now compute objectives and contraints together
                allQuantitiesIdx = unique([this.useStructsOptimization, this.useStructsConstraint]);
                this.d(allQuantitiesIdx) = arrayfun(@(struct) this.computeQuantity(dij,struct,w),allQuantitiesIdx, 'UniformOutput',false);
                this.wCache = w;
            end
            quantityOutput = this.d;
        end

        function gradOutput = getProjectedGradient(this,dij,fGrad,w)
            if ~isequal(this.wGradCache,w)
                this.wGrad(this.useStructsOptimization) = arrayfun(@(struct) this.projectGradient(dij,struct,fGrad,w), this.useStructsOptimization, 'UniformOutput',false);
                this.wGradCache = w;
            end
            gradOutput = this.wGrad;
        end

        function constJacobianOutput = getProjectedJacobian(this,dij,fJacob,w)
            if ~isequal(this.wConstJacobianCache,w)
                this.wJacob(this.useStructsConstraint) = arrayfun(@(struct) this.projectConstraintJacobian(dij,struct,fJacob,w), this.useStructsConstraint, 'UniformOutput',false);
                this.wGradCache = w;
            end
            constJacobianOutput = this.wJacob;
        end

        function initializeProperties(this,cst)
            this.cst   = cst; 
            this.d     = cell(size(cst,1),1);
            this.wGrad = cell(size(cst,1),1);
            this.wJacob = cell(size(cst,1),1);

            % tmpUseStructsOptimization = [];
            % tmpUseStructsConstraint = [];
            % for i=1:size(cst,1)
            %     for j=1:size(cst{i,6})
            %         if strcmp(cst{i,6}{j}.quantity, this.quantityName)
            %             if isa(cst{i,6}{j}, 'OmegaObjectives.matRad_VarianceObjective')
            %                 tmpUseStructsOptimization = [tmpUseStructsOptimization, i];
            %             elseif isa(cst{i,6}{j}, 'OmegaConstraints.matRad_VarianceConstraint')
            %                 tmpUseStructsConstraint = [tmpUseStructsConstraint, i];
            %             end
            %         end
            %     end
            % end
            % 
            % this.useStructsOptimization = unique(tmpUseStructsOptimization);
            % this.useStructsConstraint   = unique(tmpUseStructsConstraint);
        end

        function updateSubquantityOptimization(this,structsOptimization)
 
            for subQt=this.subQuantities'
                if isa(subQt{1},'matRad_ScalarQuantity')
                    subQt{1}.useStructsOptimization = structsOptimization;
                end
            end
        end

        function updateSubquantityConstraints(this,structsConstraint)
 
            for subQt=this.subQuantities'
                if isa(subQt{1},'matRad_ScalarQuantity')
                    subQt{1}.useStructsConstraint = structsConstraint;
                end
            end
        end

    end

    methods
        function set.useStructsOptimization(this,value)

            tmpUseStructs = unique([this.useStructsOptimization, value]);

            if ~isequal(tmpUseStructs, this.useStructsOptimization)
                this.useStructsOptimization = tmpUseStructs;
                this.updateSubquantityOptimization(this.useStructsOptimization);
            end

        end

        function set.useStructsConstraint(this,value)

            tmpUseStructs = unique([this.useStructsConstraint, value]);
            
            if ~isequal(tmpUseStructs, this.useStructsConstraint)
                this.useStructsConstraint = tmpUseStructs;
                this.updateSubquantityConstraints(this.useStructsConstraint);
            end

        end
    end
end