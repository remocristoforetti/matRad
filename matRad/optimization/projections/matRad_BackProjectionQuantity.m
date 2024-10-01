classdef matRad_BackProjectionQuantity < handle
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
        d
        wGrad
        optimizationQuantitiesIdx;
    end
    
    properties 
        dij                     %reference to matRad dij struct (to enable local changes)
        scenarios    = 1        %Scenario indices to evaluate (used for 4D & robust/stochastic optimization)
        scenarioProb = 1        %Probability associated with scenario (for stochastic optimization)
        nominalCtScenarios = 1; %nominal ct scenario (no shift, no range error) indices to evaluate (used for 4D & robust/stochastic optimization, when at least one cst structure does not have robustness)
        quantities;             % Quantities that need to be evaluated (includes subquantities)
        optimizationQuantities; % Quantities on which an objective function is defined
        structsForScalarQuantity;
    end

    
    methods
        function obj = matRad_BackProjectionQuantity()
            obj.wCache = [];
            obj.wGradCache = [];
            obj.d = [];
            obj.wGrad = [];
            obj.quantities = {};
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
        
        % function obj = computeGradientProb(obj,dij,doseGrad,vOmegaGrad,w)
        %     if ~isequal(obj.wGradCacheProb,w)
        %         obj.wGradProb = obj.projectGradientProb(dij,doseGrad,vOmegaGrad,w);
        %         obj.wGradCacheProb = w;
        %     end
        % end
        
        function d = GetResult(obj)
            d = obj.d;
        end
        
        % function [dExp,dOmegaV] = GetResultProb(obj)
        %     dExp = obj.dExp;
        %     dOmegaV = obj.dOmegaV;
        % end

        function wGrad = GetGradient(obj)
            wGrad = obj.wGrad;
        end
      
        function computeResult(obj,dij,w)
            tmpQuantitiesOutput = [];
            
            for quantityIdx=obj.optimizationQuantitiesIdx'
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
      
        function instantiateQuatities(this, optimizationQuantities, dij,cst)
            
            matRad_cfg = MatRad_Config.instance();

            if ~iscell(optimizationQuantities)
                matRad_cfg.dispError('Input quantities should be a cell array');
            end

            if ~iscolumn(optimizationQuantities)
                optimizationQuantities = optimizationQuantities';
            end
            availableQuantitiesMeta = this.getAvailableOptimizationQuantities();

            if ~all(ismember(optimizationQuantities, {availableQuantitiesMeta.quantityName}))
                matRad_cfg.dispError('Unrecognized quantity:%s',optimizationQuantities{ismember(optimizationQuantities, {availableQuantitiesMeta.quantityName})});
            end

            subQuantitiesName = {};
            for quantityIdx=1:numel(optimizationQuantities)
                currQuantityName = optimizationQuantities{quantityIdx};
                currQuantityIdx = find(ismember({availableQuantitiesMeta.quantityName}, currQuantityName));

                if ~isempty(currQuantityIdx)
                    currQuantity = availableQuantitiesMeta(currQuantityIdx).handle();
                    subQuantitiesName = [subQuantitiesName;this.getSubQuantities(currQuantity,availableQuantitiesMeta)];
                end
            end

            allQuantitiesName = [optimizationQuantities; subQuantitiesName];
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
            % [~,optimizationQuantitiesIdx] = intersect({selectedQuantitiesMeta.quantityName},optimizationQuantities);
            % this.optimizationQuantities = this.quantities(optimizationQuantitiesIdx);
        end

        
        
        % function instantiateQuatities(this, optimizationQuantities)
        % 
        %     % Get all available quantities
        %     availableQuantities = this.getAvailableOptimizationQuantities();
        % 
        %     % Find all required nested subquantities
        %     subQuantities = {};
        %     for quantityIdx=1:numel(optimizationQuantities)
        %         currQuantityName = optimizationQuantities{quantityIdx};
        %         currQuantityIdx = find(ismember({availableQuantities.quantityName}, currQuantityName));
        % 
        %         if ~isempty(currQuantityIdx)
        %             currQuantity = availableQuantities(currQuantityIdx).handle();
        %             subQuantities = [subQuantities,this.getSubQuantities(currQuantity)];
        %         end
        %     end
        % 
        %     %Convert handles to stings, easier
        %     subQuantitiesString = cellfun(@(x) func2str(x), subQuantities, 'UniformOutput',false);
        % 
        %     subQuantitiesClassName = unique(subQuantitiesString);
        % 
        %     % Make subquantities unique
        %     % tmp_subQuantities = [];
        %     % for subIdx=1:numel(subQuantities)
        %     %     repeatedSubQuantity = cellfun(@(x) isequal(subQuantities{subIdx},x), subQuantities(1:subIdx-1));
        %     %     if all(~repeatedSubQuantity)
        %     %         tmp_subQuantities = [tmp_subQuantities, subQuantities(subIdx)];
        %     %     end
        %     % end
        % 
        % 
        %     % Get the actual class handles
        %     % subQuantities = [];
        %     % for subIdx=1:numel(tmp_subQuantities)
        %     %     subQuantities = [subQuantities, {tmp_subQuantities{subIdx}()}];
        %     % end
        %     % mainQuantitiesIdx = ismember({availableQuantities.quantityName}, optimizationQuantities);
        % 
        %     % mainQuantities = [];
        %     % for qtIdx=find(mainQuantitiesIdx)
        %     %     mainQuantities = [mainQuantities,{availableQuantities(qtIdx).className}];
        %     % end
        % 
        %     mainQuantitiesClassName = {availableQuantities(ismember({availableQuantities.quantityName}, optimizationQuantities)).className}';
        % 
        % 
        %     % Get all the quantities class names
        %     allQuantitiesClassName = [mainQuantitiesClassName; subQuantitiesClassName];
        % 
        %     % Get instances
        %     allQuantitiesHandles   = cellfun(@(x) str2func(x), allQuantitiesClassName, 'UniformOutput',false);
        %     allQuantitiesInstances = cellfun(@(x) x(), allQuantitiesHandles, 'UniformOutput',false);
        % 
        %     baseQuantitiesInstances  = allQuantitiesInstances(cellfun(@(x) isempty(x.subQuantities), allQuantitiesInstances));
        %     baseQuantitiesClassNames = cellfun(@(x) class(x), baseQuantitiesInstances, 'UniformOutput',false);
        % 
        %     dependentQuantitiesInstances  = allQuantitiesInstances(cellfun(@(x) ~isempty(x.subQuantities), allQuantitiesInstances));
        %     dependentQuantitiesClassNames = cellfun(@(x) class(x),dependentQuantitiesInstances, 'UniformOutput',false);
        % 
        %     % Instantiate base quantities
        %     this.quantities = baseQuantitiesInstances;
        % 
        %     quantitiesToInstantiate = dependentQuantitiesInstances;
        %     while ~isempty(quantitiesToInstantiate)
        %         for qtIdx=1:numel(quantitiesToInstantiate)
        %             currSubQuantitiesClassNames = cellfun(@(x) func2str(x), quantitiesToInstantiate{qtIdx}.subQuantities, 'UniformOutput',false)';
        %             nonBaseSubquantities = currSubQuantitiesClassNames(~ismember(currSubQuantitiesClassNames,baseQuantitiesClassNames));
        % 
        %             if isempty(nonBaseSubquantities)
        %                 %Then I can instantiate the class with the
        %                 %subclasses and store it
        %                 currQuantityInstance = quantitiesToInstantiate{qtIdx};
        %                 [~,currBaseQuantityIdx]  = intersect(this.quantities,);
        %                 currQuantityInstance.subQuantities = %cellfun(@(x) str2func(x),currSubQuantitiesClassNames, 'UniformOutput',false);
        %                 this.quantities = [this.quantities, currQuantityInstance];
        %             end
        %         end
        % 
        %         quantitiesToInstantiateClassNames = cellfun(@(x) class(x), quantitiesToInstantiate, 'UniformOutput',false);
        %         currQuantitiesClassNames          = cellfun(@(x) class(x), this.quantities, 'UniformOutput',false);
        %         quantitiesToInstantiate = quantitiesToInstantiate(~intersect(this.quantities, quantitiesToInstantiate));
        %     end
        % 
        %     % % Store already the instances that do not have sub-quantities
        %     % baseQuantities = allQuantities(cellfun(@(x) isempty(x.subQuantities), allQuantities));
        %     % this.quantities = baseQuantities;
        %     % 
        %     % nonBaseQuantities = allQuantities(cellfun(@(x) ~isempty(x.subQuantities), allQuantities));
        %     % quantitiesToInstantiate = nonBaseQuantities;
        %     % while ~isempty(quantitiesToInstantiate)
        %     %     for qtIdx=1:numel(quantitiesToInstantiate)
        %     %         currSubquantities = quantitiesToInstantiate{qtIdx}.subQuantities;
        %     %         for subIdx=1:numel(currSubquantities)
        %     %             if ~any(cellfun(@(x) isequal(currSubquantities{subIdx}(), x), this.quantities))
        %     %                 this.quantities = [this.quantities, {currSubquantities{subIdx}()}];
        %     %             end
        %     %         end
        %     %     end
        %     %     %quantitiesToInstantiate = quantitiesToInstantiate(cellfun(@(x) ~isequal(this.quantities, x), quantitiesToInstantiate));
        %     % 
        %     % end


            
        %end
    
        function updateScenariosForQuantities(this)
            % For the time being just mirror the scenarios here, then
            % assign according to selection of the objectives
            for quantityIdx=1:numel(this.quantities)
                if isa(this.quantities{quantityIdx}, 'matRad_DistributionQuantity')
                    this.quantities{quantityIdx}.useScenarios = this.scenarios;
                end
            end
        end

        function updateStructsForQuantities(this)
            % For the time being just mirror here, then
            % assign according to selection of the objectives
            for quantityIdx=1:numel(this.quantities)
                if isa(this.quantities{quantityIdx}, 'matRad_ScalarQuantity')
                    this.quantities{quantityIdx}.useStructs = this.structsForScalarQuantity;
                end
            end
        end
    end
   
    
    %These should be abstract methods, however Octave can't parse them. As soon 
    %as Octave is able to do this, they should be made abstract again 
    % methods %(Abstract)
    %     function d = computeSingleScenario(obj,dij,scen,w)
    %         error('Function needs to be implemented');
    %     end
    % 
    %     function wGrad = projectSingleScenarioGradient(obj,dij,doseGrad,scen,w)
    %         error('Function needs to be implemented');
    %     end
    % 
    %     function [dExp,dOmegaV] = computeSingleScenarioProb(obj,dij,scen,w)
    %         %warning('');
    %     end
    % 
    %     function [dExp,dOmegaV] = projectSingleScenarioGradientProb(obj,dij,dExpGrad,dOmegaVgrad,scen,w)
    %         %warning('');
    %     end
    % end

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

                    subSubQuantities = getSubQuantities@matRad_BackProjectionQuantity(currQuantityInstance,availableQuantitiesMeta);
                    subQuantity = [subQuantity; subSubQuantities];

                end
            end

        end
    end

    methods
        function set.scenarios(this,value)
            this.scenarios = value;
            this.updateScenariosForQuantities();
        end

         function set.structsForScalarQuantity(this,value)
            this.structsForScalarQuantity = value;
            this.updateStructsForQuantities();
        end
    end
end

