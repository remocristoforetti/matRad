classdef matRad_LexicographicList < handle

    properties
        objectives;
        constraints;
    end

    properties (SetAccess=private, GetAccess=private)
        adaptedGoal=false;
    end
   

    methods
        
        function this = matRad_LexicographicList()

        end

        function addObjective(this, varargin)
        % Input can be:
        % - a LexicographicObjective object
        % The list: priority, objective, cstIdx, quantity, goal, robustness

            matRad_cfg = MatRad_Config.instance();

            if ~isa(varargin{1}, 'matRad_LexicographicObjective')
                % Create an objective from input. asume correct ordering of
                % inputs
                lexObjective = matRad_LexicographicObjective(varargin{:});
            else
                lexObjective = varargin{1};
            
            end
         
            if isempty(this.objectives) % If empty directly add the objective
                this.objectives = lexObjective;
            else
                % Add new priority to the existing ones
                currPriorities = [this.objectives.priority, lexObjective.priority];

                % Get sorting order for the objectives
                [~,sortIdx] = sort(currPriorities);

                % Add the objective
                this.objectives = [this.objectives, lexObjective];

                % Sort the order
                this.objectives = this.objectives(sortIdx);

                % If new objective is added, need to "clean" all the
                % objectives with lower priority
                if isempty(lexObjective.optimValue)
                    lowerPriorityObjectives = find([this.objectives.priority]> lexObjective.priority);

                    if ~isempty(lowerPriorityObjectives)

                        if any([this.objectives(lowerPriorityObjectives).optimValue])
                            matRad_cfg.dispWarning('Added objective has a priority higher than already optimized objectives.');
                            for i=lowerPriorityObjectives
                                this.objectives(i).resetOptimValue();
                            end
                        end
                    end
                end

            end
        end

        function addConstraint(this, varargin)
            
            matRad_cfg = MatRad_Config.instance();

            if ~isa(varargin{1}, 'matRad_LexicographicObjective')
                % Create an objective from input. asume correct ordering of
                % inputs
                if nargin == 5
                    lexObjective = matRad_LexicographicObjective([],varargin{1:3},[],varargin{4:end}); % Priority and goals not ste for constarints. Maintained if 
                else
                    lexObjective = matRad_LexicographicObjective(varargin{:});
                end
            else
                lexObjective = varargin{1};
            end
         
            if isa(lexObjective.objectiveFunction, 'DoseConstraints.matRad_DoseConstraint') || isa(lexObjective.objectiveFunction, 'OmegaConstraints.matRad_VarianceConstraint')
                this.constraints = [this.constraints, lexObjective];
            else
                matRad_cfg.dispError('Unable to assign constraint to list');
            end
        end


        function removeObjective(this, objectiveIdx)
            
            this.objectives(objectiveIdx) = [];

        end

        function finalizePriorities(this)

            % Make sure that priorities are updated correctly and are in
            % incressing order without gaps.
            currPriorities = [this.objectives.priority];
            
            [~, ~, new_idx] = unique(currPriorities);

            for i=1:numel(this.objectives)
                this.objectives(i).priority = new_idx(i);
            end
        end


        function [objectives, idxs] = getObjectivesWithPriority(this, priority)

            idxs = find([this.objectives.priority] == priority);
            objectives = [this.objectives(idxs)];
            
        end

        function adaptToFractionSize(this, numOfFractions)

            if ~this.adaptedGoal
                arrayfun(@(x) x.setGoalValue(x.objectiveFunction.adaptGoalToFraction(x.goal,numOfFractions)), this.objectives);
                this.adaptedGoal = true;
            end

        end

        function newList = createListCopy(this)

            % Create a new class instance. This could be properly done by
            % inheriting from copyable handle classes and so on, but I
            % wanna keep this simple

            newList = matRad_LexicographicList();

            % Objectives
            for i=1:numel(this.objectives)
                currObjective = this.objectives(i);

                newObjective = matRad_LexicographicObjective(currObjective.priority,currObjective.objectiveFunction,currObjective.cstIdx, currObjective.quantity,currObjective.goal, currObjective.robustness);
                newObjective.setOptimValue(currObjective.optimValue);
                
                newList.addObjective(newObjective);

            end

            % Constraints
            for i=1:numel(this.constraints)
                currObjective = this.constraints(i);
                
                if isa(currObjective.objectiveFunction, 'DoseConstraints.matRad_DoseConstraintFromObjective') || isa(currObjective.objectiveFunction, 'OmegaConstraints.matRad_OmegaConstraintFromObjective')
                    newObjective = matRad_LexicographicObjective(currObjective.priority,currObjective.objectiveFunction,currObjective.cstIdx, currObjective.quantity,currObjective.goal, currObjective.robustness);
                    newObjective.setOptimValue(currObjective.optimValue);
                else
                    newObjective = matRad_LexicographicObjective([],currObjective.objectiveFunction,currObjective.cstIdx, currObjective.quantity,[], currObjective.robustness);
                end
                newList.addConstraint(newObjective);

            end

            newList.setAdaptedGoal(this.adaptedGoal);

        end

        function structList = convertToStruct(this)
            % Converts the list to a struct for saving

            structList = struct();

            for i=1:numel(this.objectives)

                currStructList = struct();
                
                currentCstObj = this.objectives(i).objectiveFunction;

                metaClass = metaclass(currentCstObj);

                cFields = properties(currentCstObj);

                % Get rid of 'name' field
                cFields(strcmp(cFields, 'name')) = [];

                currStructObjective = struct();

                for j=1:numel(cFields)
                    currStructObjective = setfield(currStructObjective, cFields{j}, currentCstObj.(cFields{j}));
                end

                % This might be usefull to convert struct to class back
                currStructObjective.className = metaClass.Name;
                
                currStructList.objectiveFunction = currStructObjective;

                cFields = properties(this.objectives(i));
                
                for j=1:numel(cFields)
                    if ~strcmp(cFields{j}, 'objectiveFunction')
                        currStructList = setfield(currStructList, cFields{j}, this.objectives(i).(cFields{j}));
                    end
                end

                structList.objectives(i) = currStructList;
            end

            % Constaraints
            for i=1:numel(this.constraints)

                currStructList = struct();
                currentCstObj = this.constraints(i).objectiveFunction;

                metaClass = metaclass(currentCstObj);
                cFields = properties(currentCstObj);

                % Get rid of 'name' field
                cFields(strcmp(cFields, 'name')) = [];

                currStructObjective = struct();

                for j=1:numel(cFields)
                    currStructObjective = setfield(currStructObjective, cFields{j}, currentCstObj.(cFields{j}));
                end

                % This might be usefull to convert struct to class back
                currStructObjective.className = metaClass.Name;
                currStructList.objectiveFunction = currStructObjective;

                cFields = properties(this.constraints(i));
                
                for j=1:numel(cFields)
                    if ~strcmp(cFields{j}, 'objectiveFunction')
                        currStructList = setfield(currStructList, cFields{j}, this.constraints(i).(cFields{j}));
                    end
                end

                structList.constraints(i) = currStructList;
            end

            structList.adaptedGoal = this.adaptedGoal;
        end

        function loadListFromStep(this, filePath)
            
            structList = load(filePath, 'list');
            structList = structList.list;

            for i=1:numel(structList.objectives)
                currObjectiveFunction = matRad_DoseOptimizationFunction.createInstanceFromStruct(structList.objectives(i).objectiveFunction);
                objective = matRad_LexicographicObjective(structList.objectives(i).priority,currObjectiveFunction,structList.objectives(i).cstIdx,structList.objectives(i).quantity, structList.objectives(i).goal, structList.objectives(i).robustness);

                objective.setOptimValue(structList.objectives(i).optimValue);

                this.addObjective(objective);
            end

            if isfield(structList, 'constraints')
                for i=1:numel(structList.constraints)
                    currObjectiveFunction = matRad_DoseOptimizationFunction.createInstanceFromStruct(structList.constraints(i).objectiveFunction);
                    constraint = matRad_LexicographicObjective(structList.constraints(i).priority,currObjectiveFunction,structList.constraints(i).cstIdx,structList.constraints(i).quantity, structList.constraints(i).goal, structList.constraints(i).robustness);
    
                    constraint.setOptimValue(structList.constraints(i).optimValue);
    
                    this.addConstraint(constraint);
                end
            end

            this.setAdaptedGoal(structList.adaptedGoal);
        end

        function setAdaptedGoal(this, value)
            this.adaptedGoal = value;
        end

        function value = getAdaptedGoal(this)
            value = this.adaptedGoal;
        end
  
    end

    methods (Static)

        function reorderedList = reshapeList(list)
            % This functon turns back all const raints that were originaly
            % objectives into objectives again and preserves the optimized
            % values.

            % Get the objectives
            objectives = list.objectives;

            % Get constrtaints that were objectives

            constraints = list.constraints;
            for i=1:numel(list.constraints)
                if isa(list.constraints(i).objectiveFunction, 'DoseConstraints.matRad_ConstraintsFromObjective') || isa(list.constraints(i).objectiveFunction, 'OmegaConstraints.matRad_OmegaConstraintFromObjective')
                    nObjective = matRad_LexicographicObjective(); % Build a new objetive
                    nObjective.objectiveFunction = list.constraints(i).objectiveFunction.objective;
                    nObjective.cstIdx            = list.constraints(i).cstidx;
                    nObjective.priority          = list.constraints(i).priority;
                    nObjective.goal              = list.constraints(i).goal;
                    nObjective.quantity          = list.constraints(i).quantity;
                    nObjective.robustness        = list.constraints(i).robustness;

                    nObjective.setOptimValue = list.constraints(i).objectiveFunction.parameters{1};
                    % Add it to the list
                    objectives = [objectives, nObjective];

                    constraints(i) = [];
                end
            end

            % Maybe reordering not even necessary
            % Sort priorities
            [~, sIdx] = sort([objectives.priorities]);

            % Sort objectives
            objectives = [objectives(sIdx)];

            reorderedList = matRad_LexicographicList();

            for i=1:numel(objectives)
                reorderedList.addObjective(objectives(i));
            end

            for i=1:numel(constraints)
                reorderedList.addConstraint9(constraints(i));
            end
        end
    end
end