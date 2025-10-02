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
        d
        wGrad
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
    end

    
    methods
        function obj = matRad_BackProjection()

            obj.d = [];
            obj.wGrad = [];
            obj.quantities = {};
            obj.wJacob = [];

        end       
        
        function obj = compute(obj,dij,w)
            tmpQuantitiesOutput = [];
            
            allQuantitiesIdx = unique([obj.optimizationQuantitiesIdx; obj.constrainedQuantitiesIdx]);
            for quantityIdx=allQuantitiesIdx'
                quantity = obj.quantities{quantityIdx};
                tmpQuantitiesOutput.(quantity.quantityName) = quantity.getResult(dij,w);
            end
            
            obj.d = tmpQuantitiesOutput;
        end
        

        function obj = computeGradient(obj,dij,fGrad,w)
            tmpGradient = [];
            for quantityIdx=obj.optimizationQuantitiesIdx'
                quantity = obj.quantities{quantityIdx};
                tmpGradient.(quantity.quantityName) = quantity.getProjectedGradient(dij,fGrad.(quantity.quantityName),w);
            end
            obj.wGrad = tmpGradient;
        end

        function obj = computeConstraintJacobian(obj,dij,fJacob, w)
            % if ~isequal(obj.wConstJacobianCache,w)
            %     obj.projectConstraintJacobian(dij,fJacob,w);
            %     obj.wConstJacobianCache = w;
            % end
            tmpGradient = [];
            for quantityIdx=obj.constrainedQuantitiesIdx'
                quantity = obj.quantities{quantityIdx};
                tmpGradient.(quantity.quantityName) = quantity.getProjectedJacobian(dij,fJacob.(quantity.quantityName),w);
            end
            obj.wJacob = tmpGradient;
        end
        
        function d = GetResult(obj)
            d = obj.d;
        end

        function wGrad = GetGradient(obj)
            wGrad = obj.wGrad;
        end
      
        function instantiateQuatities(this, optimizationQuantities, constraintQuantities, dij)
            
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
                constraintQuantities = constraintQuantities(:);
            end
 
            availableQuantitiesMeta = this.getAvailableOptimizationQuantities();

            % Check for unrecognized provided quantities
            if ~all(ismember(optimizationQuantities, {availableQuantitiesMeta.quantityName}))
                matRad_cfg.dispError('Unrecognized quantity:%s',optimizationQuantities{~ismember(optimizationQuantities, {availableQuantitiesMeta.quantityName})});
            end

            if ~all(ismember(constraintQuantities, {availableQuantitiesMeta.quantityName}))
                matRad_cfg.dispError('Unrecognized quantity:%s',constraintQuantities{~ismember(constraintQuantities, {availableQuantitiesMeta.quantityName})});
            end

            % Check for subquantiites of the requested optimized/constrained quantities
            allQuantities = unique([optimizationQuantities; constraintQuantities]);
            subQuantitiesName = {};
            for quantityIdx=1:numel(allQuantities)
                currQuantityName = allQuantities{quantityIdx};
                currQuantityIdx  = find(ismember({availableQuantitiesMeta.quantityName}, currQuantityName));

                if ~isempty(currQuantityIdx)
                    currQuantity      = availableQuantitiesMeta(currQuantityIdx).handle();
                    subQuantitiesName = [subQuantitiesName;this.getSubQuantities(currQuantity,availableQuantitiesMeta)];
                end
            end

            allQuantitiesName = [allQuantities; subQuantitiesName];
            allQuantitiesName = unique(allQuantitiesName);
            
            selectedQuantitiesMeta = availableQuantitiesMeta(ismember({availableQuantitiesMeta.quantityName}, allQuantitiesName));
            
            %Instantiate the quantities
            this.quantities = cellfun(@(x) x(), {selectedQuantitiesMeta.handle}, 'UniformOutput',false)';
            distributionQuantities = cellfun(@(x) isa(x, 'matRad_DistributionQuantity'), this.quantities);
            
            if  any(distributionQuantities)
                cellfun(@(x) x.initializeProperties(dij), this.quantities(distributionQuantities));
            end
            
            % Link quantities and subquantities
            for quantityIdx=1:numel(this.quantities)
                requiredSubquantitiesName    = this.quantities{quantityIdx}.requiredSubquantities;
                [~,requiredSubquantitiesIdx] = intersect({selectedQuantitiesMeta.quantityName},requiredSubquantitiesName);
                
                if ~isempty(requiredSubquantitiesIdx)
                    this.quantities{quantityIdx}.subQuantities = this.quantities(requiredSubquantitiesIdx);
                end
            end

            % Select the quantities that are explicitly optimized/constrained.
            this.optimizationQuantities = optimizationQuantities';
            [~,this.optimizationQuantitiesIdx] = intersect({selectedQuantitiesMeta.quantityName},optimizationQuantities);
            this.constrainedQuantities = constraintQuantities';
            [~, this.constrainedQuantitiesIdx] = intersect({selectedQuantitiesMeta.quantityName},constraintQuantities);
            
            % for i=1:size(cst,1)
            %     for j=1:size(cst{i,6},2)
            %         obj = cst{i,6}{j};
            %         qtIdx = find(arrayfun(@(qtMeta) strcmp(qtMeta.quantityName, obj.quantity), selectedQuantitiesMeta));
            %         if isa(this.quantities{qtIdx}, 'matRad_ScalarQuantity')
            %             if isa(obj, 'DoseObjectives.matRad_DoseObjective') || isa(obj, 'OmegaObjectives.matRad_OmegaObjective')
            %                 this.quantities{qtIdx}.useStructsOptimization = [this.quantities{qtIdx}.useStructsOptimization,i];
            %             elseif isa(obj, 'DoseConstraints.matRad_DoseConstraint') || isa(obj, 'OmegaConstraints.matRad_VarianceConstraint')
            %                 this.quantities{qtIdx}.useStructsConstraint = [this.quantities{qtIdx}.useStructsConstraint,i];
            %             end
            %         end
            % 
            %         if any(cellfun(@(qt) isa(qt, 'matRad_ScalarQuantity'), this.quantities{qtIdx}.subQuantities))
            %             subQtIdx = find(cellfun(@(qt) isa(qt, 'matRad_ScalarQuantity'), this.quantities{qtIdx}.subQuantities));
            %             subQtIdx = subQtIdx(:)';
            %             for k = subQtIdx
            %                 if isa(obj, 'DoseObjectives.matRad_DoseObjective') || isa(obj, 'OmegaObjectives.matRad_OmegaObjective')
            %                     this.quantities{qtIdx}.subQuantities{k}.useStructsOptimization = [this.quantities{qtIdx}.subQuantities{k}.useStructsOptimization,i];
            %                 elseif isa(obj, 'DoseConstraints.matRad_DoseConstraint') || isa(obj, 'OmegaConstraints.matRad_VarianceConstraint')
            %                     this.quantities{qtIdx}.subQuantities{k}.useStructsConstraint = [this.quantities{qtIdx}.subQuantities{k}.useStructsConstraint,i];
            %                 end
            % 
            %             end
            %         end
            % 
            %     end
            % end
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

        function w = initializeWeights(this,cst,dij, wInit)

            matRad_cfg = MatRad_Config.instance();


            if exist('wInit', 'var') && ~isempty(wInit)
                w = wInit;
                return;
            end

            % Get all Taget indexes
            allTargetIdx = find(strcmp(cst(:,3), 'TARGET'));

            % Only keep the last opne for simplicity
            for tmpTidx=allTargetIdx'
                if ~isempty(cst{tmpTidx,6})
                    targetIdx = tmpTidx;
                end
            end
                
            % Get quantities for all structures.
            allTargetQuantities = cellfun(@(x) x.quantity, cst{targetIdx,6}, 'UniformOutput', false);

            % I want to initialize the weights on the primary
            % quantities only if they are available. Example:
            % PhysicalDose + LET objectives on the target, I only want
            % to consider the physical dose for the weight
            % initialization
            primaryQuantities = {'physicalDose', 'constantRBExDose','RBExDose','BED','effect'};

            targetQtForInitialization = intersect(allTargetQuantities, primaryQuantities);
                
            if isempty(targetQtForInitialization)
                % If none of the primary quantities is set, just take the first one
                targetQtForInitialization = allTargetQuantities(1);
            end

            matRad_cfg.dispInfo(sprintf('%s selected as quantity for the weight initialization\n', targetQtForInitialization{:}));

            % Initialize the weights
            wOnes = ones(dij.totalNumOfBixels,1);

            % Get prescription
            % Select all dose objectives/constraints
            targetObjIdx = find(cellfun(@(x) (isa(x, 'DoseObjectives.matRad_DoseObjective') || ...
                                         isa(x, 'DoseConstraints.matRad_DoseConstraint')) ...
                                        && any(strcmp(x.quantity, targetQtForInitialization)),...    % And be a selected quantity for initialization
                                        cst{targetIdx,6}));
                            
            % Only catch the "effect" case, as it's the only case with
            % a different strategy

            % targetQtDorOpt, not in the same order as qts in cst
            % becasue of intersect function
            bixelWeight = [];
            for objIdx=targetObjIdx
                
                currTargetQuantity = cst{targetIdx,6}{objIdx}.quantity;
                
                if strcmp(currTargetQuantity, 'effect')

                    % This is to reproduce old branch
                    targetObjective = matRad_Effect.setBiologicalDosePrescriptions(cst{targetIdx,6}{targetObjIdx}, cst{targetIdx,5}.alphaX, cst{targetIdx,5}.betaX);
                    d_pres_effect = targetObjective.getDoseParameters;

                    d_pres_effect = d_pres_effect(~isinf(d_pres_effect));
                    d_pres_effect = mean(d_pres_effect);

                    % Set biological dose prescription
                    % Get alpha/beta instances. No photon effect for now.
                    % Easy to add later
                    try

                        alphaQuantity    = this.getQuantityNamed('AlphaDose');
                        sqrtBetaQuantity = this.getQuantityNamed('SqrtBetaDose');
                    catch
                        try 
                            % For probabilistic planning
                            alphaQuantity    = this.getQuantityNamed('AlphaDoseExp');
                            sqrtBetaQuantity = this.getQuantityNamed('SqrtBetaDoseExp');
                        catch

                        end
                    end

                    aTmp = alphaQuantity.computeQuantity(dij,1,wOnes);
                    bTmp = sqrtBetaQuantity.computeQuantity(dij,1,wOnes);

                    % This cst is already on the grid
                    V = cst{targetIdx,4}{1};
                    p = sum(aTmp(V))/sum(bTmp(V).^2);
                    q = -(d_pres_effect*length(V))/sum(bTmp(V).^2);

                    bixelWeight = -p/2 + sqrt((p^2)/4 - q);

                else

                    % The bixelWeight is just: TargetPrescription/mean(qt_wOnes)
                    
                    % get the average out of all the parametrs
                    d_pres_qt = cst{targetIdx,6}{objIdx}.getDoseParameters();
                    d_pres_qt = mean(d_pres_qt(~isinf(d_pres_qt)));
                    
                    % Get the taget voxels
                    V = unique(cat(1,cst{targetIdx,4}{:})); % Getting all CT scenarios here, but d computed only on the first scenario later

                    % Instantiate the quantities. targetQtForInitilaization
                    % is now unique cell, but could change in the future
                    qtInstances = this.getQuantityNamed(currTargetQuantity);
                    
                    % Compute the quantity with the wOnes weigths
                    % d_wOnes = cellfun(@(qT) qT.computeQuantity(dij, 1,wOnes), qtInstances, 'UniformOutput',false);
                    d_wOnes = qtInstances.computeQuantity(dij, 1,wOnes);

                    % If d is a distribution, get the voxels, if its a
                    % scalar, just get the value
                    if isa(qtInstances, 'matRad_DistributionQuantity')
                        d_wOnes = mean(d_wOnes(V));
                    end
                    
                    % Get bixels weight
                    bixelWeight = [bixelWeight, d_pres_qt/d_wOnes];

                end
            end
            if isempty(bixelWeight)
                bixelWeight = 1;
            end
            
            w = mean(bixelWeight,2) * wOnes;

        end

        function qtInstance = getQuantityNamed(this, qtName)
            
            % Get the quantities that are optimized
            availableQuantities = this.quantities;
            availableQuantitiesNames = cellfun(@(qt) qt.quantityName, availableQuantities, 'UniformOutput', false);

            qtIdx = find(strcmp(availableQuantitiesNames, qtName));

            if isempty(qtIdx)
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError(sprintf('Quantity called: %s not instantiated', qtName));
            end

            qtInstance = availableQuantities{qtIdx};

        end



    end
   

    methods (Static)

        % function optiFunc = setBiologicalDosePrescriptions(optiFunc,alphaX,betaX)
        %     %Does nothing in a usual normal setting but return the original
        %     %optiFunc
        % end

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

        function quantityInstance = getQuantityInstanceFromName(qtName)
            
            availableQuantitiesMeta = matRad_BackProjection.getAvailableOptimizationQuantities();

            [~,idx] = intersect({availableQuantitiesMeta(:).quantityName}, qtName);

            if ~isempty(idx)
                quantityInstance = availableQuantitiesMeta(idx).handle();
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError(sprintf('%s not recognized as an available quantity', qtName));
                quantityInstance = [];
            end
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

        function [optQuantities,constraintQuantities] = getOptimizationConstraintQuantitiesFromCst(cst, defaultQuantity)
            
            matRad_cfg = MatRad_Config.instance();
           

            tmpOptimizationQt = [];
            tmpConstraintQt = [];
            
            for i=1:size(cst,1)
                for j=1:numel(cst{i,6})
                    
                   if isempty(cst{i,6}{j}.quantity)
                        matRad_cfg.dispWarning(sprintf('No quantity detected for objecive #%d of structure:%s, using default quantity:%s', j, cst{i,2}, defaultQuantity));
                        cst{i,6}{j}.quantity = defaultQuantity;
                   end

                    if isa(cst{i,6}{j}, 'DoseObjectives.matRad_DoseObjective')
                     
                        tmpOptimizationQt = [tmpOptimizationQt, {cst{i,6}{j}.quantity}];
                    
                    elseif isa(cst{i,6}{j}, 'DoseConstraints.matRad_DoseConstraint')

                        tmpConstraintQt = [tmpConstraintQt, {cst{i,6}{j}.quantity}];
                    end
                end
            end
            
            optQuantities = unique(tmpOptimizationQt);

            if isempty(optQuantities)
                optQuantities = {};
            end            
            
            constraintQuantities = unique(tmpConstraintQt);
            if isempty(constraintQuantities)
                constraintQuantities = {};
            end
        end

    %     function G = getOptimizationQuantityGraph(verbosity)
    % 
    %         matRad_cfg = MatRad_Config.instance();
    % 
    %         if ~exist('verbosity', 'var')
    %             verbosity = false;
    %         end
    % 
    %         availableQuantitiesMeta = matRad_BackProjection.getAvailableOptimizationQuantities();
    % 
    %         nodeNames = {availableQuantitiesMeta.quantityName};
    %         A = [];
    %         for curQtName=nodeNames
    %             curQt = matRad_BackProjection.getQuantityInstanceFromName(curQtName{1});
    %             curSubQuantities = curQt.requiredSubquantities;
    %             if ~isempty(curSubQuantities)
    %                 idxs = cellfun(@(x) find(strcmp(nodeNames, x)), curSubQuantities);
    %             else
    %                 idxs = [];
    %             end
    %             currConnections = zeros(size(nodeNames));
    %             currConnections(idxs) = 1;
    % 
    %             A = [A; currConnections];
    %         end
    % 
    %         G = digraph(A, nodeNames');
    %         % G.Nodes.Properties.RowNames = nodeNames;
    % 
    %         if verbosity
    %             f = figure;
    %             f.Position(1:2) = f.Position(1:2) - 500;
    %             plot(G);
    %         end
    %     end
    
    end

    methods
        function set.scenarios(this,value)
            this.scenarios = value;
            this.updateScenariosForQuantities();
        end

    end
end

