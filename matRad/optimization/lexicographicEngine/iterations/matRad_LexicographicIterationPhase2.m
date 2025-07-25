classdef matRad_LexicographicIterationPhase2 < matRad_LexicographicIteration

    properties (Constant)
        name = 'phase2';
    end

    properties
        previousPriority = 0; % Initialize this to 0
    end

    methods

        function this = matRad_LexicographicIterationPhase2(pln)

            this@matRad_LexicographicIteration(pln);

        end

        function [objective, priority] = scrollList(this,list)

            % This function needs to decide which objective in the list
            % needs to be optimized now. In Phase 2: need to keep track
            % internaly of what was optimized before.
            
            % In this phase, list can have missing priorities, those that 
            % were not achievable in phase 1 will be constraints at this
            % point.
            prioritiesInList = [0,list.objectives.priority]; % Include 0 for indexing

            % Move to next priority. This just finds the index of the
            % previous poriority in the list of priorities to optimize and
            % picks the next one.
            currPriority = prioritiesInList(find(prioritiesInList == this.previousPriority,1,'last')+1);
            
            % Find the objective(s) with current priority
            
            objective = list.getObjectivesWithPriority(currPriority);
            priority = currPriority;

            % Update the previous priority.
            this.previousPriority = currPriority;
        end

        function nIterations = getNumOfIterations(this,list)

            % All the objectives need to be optimized
            % nIterations = numel(unique([list.objectives.priority]));            
            nIterations = sum(unique([list.objectives.priority]) > this.previousPriority); % Exclude those with lower priority if it's initilized manually

           
        end

        function cst = getCst(this, list, originalCst, objective, currPriority)

            cst = originalCst;
            
            cst(:,6) = cell(size(cst,1),1);

            % Add the current objective
            for i =1:numel(objective)
                cstObj            = objective(i).objectiveFunction;
                cstObj.quantity   = objective(i).quantity;
                cstObj.robustness = objective(i).robustness;

                cst{objective(i).cstIdx, 6} = [cst{objective(i).cstIdx, 6}, {cstObj}];
                
            end

            roundObjective = @(x,n) round(x*10^n)*10^(-n);

            % Turn everything else in constraint. Here have to decide what
            % bound to apply
            for oIdx=1:numel(list.objectives) % loop over all the objectives
                currObjective = list.objectives(oIdx);
                
                if currObjective.priority ~= currPriority % skip current objective
        
                    if ~isempty(currObjective.optimValue) % If it's empty we have a problem
                       maxObj = roundObjective(currObjective.optimValue*this.slack, 8);
                    end

                    cstObj            = currObjective.objectiveFunction.turnIntoLexicographicConstraint(maxObj);
                    cstObj.quantity   = currObjective.quantity;
                    cstObj.robustness = currObjective.robustness;

                    cst{currObjective.cstIdx, 6} = [cst{currObjective.cstIdx, 6}, {cstObj}];
                    
                end
            end

            % Add constraints
            for oIdx=1:numel(list.constraints) % loop over all the objectives
                currConstraint = list.constraints(oIdx);
                
                cstObj            = currConstraint.objectiveFunction;
                cstObj.quantity   = currConstraint.quantity;
                cstObj.robustness = currConstraint.robustness;

                cst{currConstraint.cstIdx, 6} = [cst{currConstraint.cstIdx, 6}, {cstObj}];
            end

        end

        function updatedList = updateList(this, list, values, problemMeta)
            
            % Get the objectives with the current priority
            [objectives, idxs] = list.getObjectivesWithPriority(problemMeta.stepPriority);

            for i=1:numel(objectives)
               objectives(i).setOptimValue(values(i));
            end
            
            updatedList = list;
            updatedList.objectives(idxs) = objectives;
        end

        
        
   end
end