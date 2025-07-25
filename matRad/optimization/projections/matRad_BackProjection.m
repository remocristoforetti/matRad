classdef matRad_BackProjection < handle
% matRad_BackProjection superclass for all backprojection algorithms 
% used within matRad optimzation processes
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    properties (SetAccess = protected)
        wCache
        wGradCache  %different cache for optimal performance (if multiple evaluations of objective but not gradient are required)
        %wConstraintCache
        wConstJacobianCache
        d
        wGrad
        %c
        wJacob
        optimizationQuantitiesIdx;
        constrainedQuantitiesIdx;
    end
    
    properties 
        dij                     %reference to matRad dij struct (to enable local changes)
        scenarios    = 1        %Scenario indices to evaluate (used for 4D & robust/stochastic optimization)
        scenarioProb = 1        %Probability associated with scenario (for stochastic optimization)
        nominalCtScenarios = 1; %nominal ct scenario (no shift, no range error) indices to evaluate (used for 4D & robust/stochastic optimization, when at least one cst structure does not have robustness)
        quantities;             % Quantities that need to be evaluated (includes subquantities)
        optimizationQuantities; % Quantities on which an objective function is defined
        constrainedQuantities;
        %structsForScalarQuantity;
    end

    
    methods
        function obj = matRad_BackProjection()
            obj.wCache = [];
            obj.wGradCache = [];
            obj.d = [];
            obj.wGrad = [];
            obj.quantities = {};
            obj.wJacob = [];
            obj.wConstJacobianCache = [];
        end       
        
        function obj = compute(obj,dij,w)
            if ~isequal(obj.wCache,w)
                obj.computeResult(dij,w);
                obj.wCache = w;
            end
        end
        
        function obj = computeGradient(obj,dij,fGrad,w)
            if ~isequal(obj.wGradCache,w)
                obj.projectGradient(dij,fGrad,w);
                obj.wGradCache = w;
            end
        end
        
        function obj = computeConstraintJacobian(obj,dij,fJacob, w)
            if ~isequal(obj.wConstJacobianCache,w)
                obj.projectConstraintJacobian(dij,fJacob,w);
                obj.wConstJacobianCache = w;
            end
        end
        
        function d = GetResult(obj)
            d = obj.d;
        end
        
        function wGrad = GetGradient(obj)
            wGrad = obj.wGrad;
        end
      
        function computeResult(obj,dij,w)
            tmpQuantitiesOutput = [];
            
            allQuantitiesIdx = unique([obj.optimizationQuantitiesIdx; obj.constrainedQuantitiesIdx]);
            for quantityIdx=allQuantitiesIdx'
                quantity = obj.quantities{quantityIdx};
                tmpQuantitiesOutput.(quantity.quantityName) = quantity.getResult(dij,w);
            end
            obj.d = tmpQuantitiesOutput;

        end

        function projectGradient(obj,dij,fGrad,w)
            tmpGradient = [];
            for quantityIdx=obj.optimizationQuantitiesIdx'
                quantity = obj.quantities{quantityIdx};
                tmpGradient.(quantity.quantityName) = quantity.getProjectedGradient(dij,fGrad.(quantity.quantityName),w);
            end
            obj.wGrad = tmpGradient;
        end
     
        function projectConstraintJacobian(obj,dij,fJacob,w)

            tmpGradient = [];
            for quantityIdx=obj.constrainedQuantitiesIdx'
                quantity = obj.quantities{quantityIdx};
                tmpGradient.(quantity.quantityName) = quantity.getProjectedJacobian(dij,fJacob.(quantity.quantityName),w);
            end
            obj.wJacob = tmpGradient;
        end
      
        
        function instantiateQuatities(this, optimizationQuantities, constraintQuantities, dij,cst)
            
            matRad_cfg = MatRad_Config.instance();

            if ~iscell(optimizationQuantities)
                matRad_cfg.dispError('Input quantities should be a cell array');
            end

            if ~iscolumn(optimizationQuantities)
                optimizationQuantities = optimizationQuantities';
            end
  
            if ~iscell(constraintQuantities)
                matRad_cfg.dispError('Input quantities should be a cell array');
            end

            if ~iscolumn(constraintQuantities)
                constraintQuantities = constraintQuantities';
            end
 
            availableQuantitiesMeta = this.getAvailableOptimizationQuantities();

            if ~all(ismember(optimizationQuantities, {availableQuantitiesMeta.quantityName}))
                matRad_cfg.dispError('Unrecognized quantity:%s',optimizationQuantities{~ismember(optimizationQuantities, {availableQuantitiesMeta.quantityName})});
            end

            if ~all(ismember(constraintQuantities, {availableQuantitiesMeta.quantityName}))
                matRad_cfg.dispError('Unrecognized quantity:%s',constraintQuantities{~ismember(constraintQuantities, {availableQuantitiesMeta.quantityName})});
            end

            allOptConstrainedQuantities = unique([optimizationQuantities; constraintQuantities]);
            subQuantitiesName = {};
            for quantityIdx=1:numel(allOptConstrainedQuantities)
                currQuantityName = allOptConstrainedQuantities{quantityIdx};
                currQuantityIdx = find(ismember({availableQuantitiesMeta.quantityName}, currQuantityName));

                if ~isempty(currQuantityIdx)
                    currQuantity = availableQuantitiesMeta(currQuantityIdx).handle();
                    subQuantitiesName = [subQuantitiesName;this.getSubQuantities(currQuantity,availableQuantitiesMeta)];
                end
            end

            allQuantitiesName = [allOptConstrainedQuantities; subQuantitiesName];
            allQuantitiesName = unique(allQuantitiesName);
            
            selectedQuantitiesMeta = availableQuantitiesMeta(ismember({availableQuantitiesMeta.quantityName}, allQuantitiesName));
            %Instantiate the quantities
           
            this.quantities = cellfun(@(x) x(), {selectedQuantitiesMeta.handle}, 'UniformOutput',false)';
            
            distributionQuantities = cellfun(@(x) isa(x, 'matRad_DistributionQuantity'), this.quantities);
            

            if  any(distributionQuantities)
                cellfun(@(x) x.initializeProperties(dij), this.quantities(distributionQuantities));
            end

            if any(~distributionQuantities)
                cellfun(@(x) x.initializeProperties(cst), this.quantities(~distributionQuantities));
            end
            
            for quantityIdx=1:numel(this.quantities)
                requiredSubquantitiesName = this.quantities{quantityIdx}.requiredSubquantities;
                [~,requiredSubquantitiesIdx] = intersect({selectedQuantitiesMeta.quantityName},requiredSubquantitiesName);
                
                if ~isempty(requiredSubquantitiesIdx)
                    this.quantities{quantityIdx}.subQuantities = this.quantities(requiredSubquantitiesIdx);
                end
            end

            this.optimizationQuantities = optimizationQuantities';
            [~,this.optimizationQuantitiesIdx] = intersect({selectedQuantitiesMeta.quantityName},optimizationQuantities);
            this.constrainedQuantities = constraintQuantities';
            [~, this.constrainedQuantitiesIdx] = intersect({selectedQuantitiesMeta.quantityName},constraintQuantities);
            % [~,optimizationQuantitiesIdx] = intersect({selectedQuantitiesMeta.quantityName},optimizationQuantities);
            % this.optimizationQuantities = this.quantities(optimizationQuantitiesIdx);
            
            for i=1:size(cst,1)
                for j=1:size(cst{i,6},2)
                    obj = cst{i,6}{j};
                    qtIdx = find(arrayfun(@(qtMeta) strcmp(qtMeta.quantityName, obj.quantity), selectedQuantitiesMeta));
                    if isa(this.quantities{qtIdx}, 'matRad_ScalarQuantity')
                        if isa(obj, 'DoseObjectives.matRad_DoseObjective') || isa(obj, 'OmegaObjectives.matRad_OmegaObjective')
                            this.quantities{qtIdx}.useStructsOptimization = [this.quantities{qtIdx}.useStructsOptimization,i];
                        elseif isa(obj, 'DoseConstraints.matRad_DoseConstraint') || isa(obj, 'OmegaConstraints.matRad_VarianceConstraint')
                            this.quantities{qtIdx}.useStructsConstraint = [this.quantities{qtIdx}.useStructsConstraint,i];
                        end
                    end
                end
            end
        end

    
        function updateScenariosForQuantities(this)
            % For the time being just mirror the scenarios here, then
            % assign according to selection of the objectives
            for quantityIdx=1:numel(this.quantities)
                if isa(this.quantities{quantityIdx}, 'matRad_DistributionQuantity')
                    this.quantities{quantityIdx}.useScenarios = this.scenarios;
                end
            end
        end

    end


    methods (Static)
        function optiFunc = setBiologicalDosePrescriptions(optiFunc,alphaX,betaX)
            %Does nothing in a usual normal setting but return the original
            %optiFunc
        end

        function quantityInfos = getAvailableOptimizationQuantities()

            matRad_cfg = MatRad_Config.instance();

            mainFolder        = fullfile(matRad_cfg.matRadSrcRoot,'optimization', 'optimizationQuantities');
            userDefinedFolder = fullfile(matRad_cfg.primaryUserFolder, 'optimizationQuantities');

            if ~exist(userDefinedFolder,"dir")
                folders = {mainFolder};
            else
                folders = {mainFolder,userDefinedFolder};
            end
            
            availableQuantitiesClassList = matRad_findSubclasses('matRad_OptimizationQuantity', 'folders', folders , 'includeSubfolders',true);

            quantityInfos = matRad_identifyClassesByConstantProperties(availableQuantitiesClassList,'quantityName');

        end

        function subQuantity = getSubQuantities(quantity, availableQuantitiesMeta)

            
            subQuantity = quantity.requiredSubquantities';
            currLevelQuantities  = subQuantity;
            nQuantities  = numel(currLevelQuantities);

            for subIdx=1:nQuantities
                currQuantityName = currLevelQuantities{subIdx};
                [~,classIdx] = intersect({availableQuantitiesMeta.quantityName},currQuantityName);
                currQuantityInstance = availableQuantitiesMeta(classIdx).handle();
                
                if isa(currQuantityInstance, 'matRad_OptimizationQuantity')

                    subSubQuantities = getSubQuantities@matRad_BackProjection(currQuantityInstance,availableQuantitiesMeta);
                    subQuantity = [subQuantity; subSubQuantities];

                end
            end

        end

        function [optQuantities,constraintQuantities] = getOptimizationConstraintQuantitiesFromCst(cst)
            
            useStructsForOmega = [];
            useStructsForConstraintOmega = [];
            omegaQuantity = [];
            quantitiesFromCst = [];
            constraintQuantities = {};
            for i=1:size(cst,1)
                for j=1:numel(cst{i,6})
                    if isa(cst{i,6}{j}, 'DoseObjectives.matRad_DoseObjective') || isa(cst{i,6}{j}, 'OmegaObjectives.matRad_OmegaObjective')
                        if isa(cst{i,6}{j}, 'OmegaObjectives.matRad_OmegaObjective') || any(strcmp(cst{i,6}{j}.quantity, {'MeanAverageEffect', 'MeanEffect', 'meanLETd', 'meanPhysicalDose'}))
                            omegaQuantity = cst{i,6}{j}.quantity;
                            useStructsForOmega = [useStructsForOmega,i];
                        elseif isa(cst{i,6}{j}, 'DoseObjectives.matRad_DoseObjective') && isempty(cst{i,6}{j}.quantity)
                            cst{i,6}{j}.quantity = pln.propOpt.quantityOpt;
                        end
                        quantitiesFromCst = [quantitiesFromCst, {cst{i,6}{j}.quantity}];
                    elseif isa(cst{i,6}{j}, 'DoseConstraints.matRad_DoseConstraint') || isa(cst{i,6}{j}, 'OmegaConstraints.matRad_VarianceConstraint')
                        constraintQuantities = [constraintQuantities, {cst{i,6}{j}.quantity}];
                        if isa(cst{i,6}{j}, 'OmegaConstraint.matRad_VarianceConstraint')
                            useStructsForOmega = [useStructsForOmega,i];
                        end
                    end
                end
            end
            
            quantitiesFromCst = unique(quantitiesFromCst);
            optQuantities = [quantitiesFromCst, {omegaQuantity}];
            optQuantities(cellfun(@isempty,optQuantities)) = [];
            optQuantities = unique(optQuantities);
            
            constraintQuantities(cellfun(@isempty,constraintQuantities)) = [];
            constraintQuantities = unique(constraintQuantities);
        end

        function w = initializeWeights(cst,dij, wInit)

            matRad_cfg = MatRad_Config.instance();

            
            if exist('wInit', 'var') && ~isempty(wInit)
                w = wInit;
            else
                % Get the otimization quantities form cst
                %optQuantities = matRad_BackProjection.getOptimizationConstraintQuantitiesFromCst(cst);

                % Get available quantities
                availableQuantitiesMeta = matRad_BackProjection.getAvailableOptimizationQuantities();

                % Get optimization quantities in the target

                % Get all Taget indexes
                allTargetIdx = find(strcmp(cst(:,3), 'TARGET'));

                % Only keep the last opne for simplicity
                for tmpTidx=allTargetIdx'
                    if ~isempty(cst{tmpTidx,6})
                        targetIdx = tmpTidx;
                    end
                end
                
                % Get quantities for all structures
                allTargetQuantities = unique(cellfun(@(x) x.quantity, cst{targetIdx,6}, 'UniformOutput', false));

                % Filter out one of the main quantities
                if numel(allTargetQuantities)>1

                    if ~isempty(intersect(allTargetQuantities, {'physicalDose', 'physicalDoseExp','constRBE', 'effect', 'RBExDose', 'BED'}))
                        % Get the main quantity giving priority to the
                        % first one if multiple are found
                        % targetQuantity = allTargetQuantities{find(cellfun(@(qt) strcmp(qt, {'physicalDose', 'constRBE', 'effect', 'RBExDose', 'BED'}), allTargetQuantities), 1, 'last'))};
                        tmpQt = intersect(allTargetQuantities, {'physicalDose', 'physicalDoseExp','constRBE', 'effect', 'RBExDose', 'BED'});
                        targetQuantity = tmpQt{1};
                    else
                        % If no special quantity defined, take the first
                        % one
                        targetQuantity = allTargetQuantities{1};
                    end
                else
                        targetQuantity = allTargetQuantities{1};
                end

                matRad_cfg.dispInfo(sprintf('%s selected as quantity for the weight initialization\n', targetQuantity));

                % Get the class handle
                
                if ~isempty(availableQuantitiesMeta(strcmp({availableQuantitiesMeta.quantityName}, targetQuantity)))
                    qtTargetHandle = availableQuantitiesMeta(strcmp({availableQuantitiesMeta.quantityName}, targetQuantity)).handle;
                else
                    matRad_cfg.dispError(sprintf('Unrecognizerd quantity %s',targetQuantity));
                end
                
                wOnes = ones(dij.totalNumOfBixels,1);

                % Compute quantity
                d = qtTargetHandle().computeQuantity(dij,1,wOnes);

                % Select dose objectives/constraints
                targetObjIdx = cellfun(@(x) isa(x, 'DoseObjectives.matRad_DoseObjective') || isa(x, 'DoseConstraints.matRad_DoseConstraint'), cst{targetIdx,6});
                
                % Get prescriptions for all objectives for the target
                d_pres = cellfun(@(x) x.getDoseParameters(), cst{targetIdx,6}(targetObjIdx), 'UniformOutput',false);
                d_pres = cellfun(@(x) x(~isinf(x)), d_pres, 'UniformOutput',false);
                %d_pres = mean([d_pres{:}]);
                % Get maximum. This might be a bit week, there is no
                % guarantee that the reference here is even linked to the
                % objective setting the quantity. This can lead to errors,
                % but should give similar behaviour wrt current
                % implementation.
                % d_pres = max([d_pres{:}]);
                d_pres = mean([d_pres{:}]);

                % Normalize teh weights to get the value for now. Then idk
                % how to handle this further avoiding the switch cases
                w = (d_pres/mean(d(cst{targetIdx,4}{1}))) * wOnes;

                % if strcmp(targetQuantity, 'constRBEDose')
                %     w = w./dij.RBE;
                % end
            end

        end
    end

    methods
        function set.scenarios(this,value)
            this.scenarios = value;
            this.updateScenariosForQuantities();
        end
    
    end
end