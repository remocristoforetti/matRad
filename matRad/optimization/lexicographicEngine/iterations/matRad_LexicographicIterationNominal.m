classdef matRad_LexicographicIterationNominal < matRad_LexicographicIteration

    properties (Constant)
        name = 'nominal';
    end

    methods

        function this = matRad_LexicographicIterationNominal(pln)

            this@matRad_LexicographicIteration(pln);

        end

        function [objective, priority] = scrollList(this,list)

           % Just scroll through the list and get the next objective to
           % optimize
            
            nextObjective = [];
            i=1;
            
            while isempty(nextObjective) % Keep going untill you find something

                if isempty(list.objectives(i).optimValue) % If this obj has not been optimized

                    priority = list.objectives(i).priority; % Get priority

                    nextObjective = list.getObjectivesWithPriority(priority); % These become the next objectives to be optimized
                else
                    i = i+1; % Move to next objective
                end

                % Check for list overflow
                if i>numel(list.objectives)
                    matRad_cfg = MatRad_Config.instance();

                    matRad_cfg.dispWarning('All objectives in the list have been optimized')

                    nextObjective = NaN;
                end

            end

            objective = nextObjective;

        end

        function nIterations = getNumOfIterations(this,list)

            % get all the list steps
            nPriorities = numel(unique([list.objectives.priority]));

            % Need to subtract the ones already optimized
            
            nIterations = -1;
            i=0;
            while nIterations == -1

               objectives = list.getObjectivesWithPriority(i+1);

               if ~isempty(objectives)
                   if any(isempty([objectives.optimValue])) % If values have not been optimized, anything below that priority will be optimized
                        nIterations = nPriorities - i;
                   else
                       i= i+1;

                   end
               else
                   % Skipping one priority?
                    i = i +1;
               end

               if i>= nPriorities
                   matRad_cfg = MatRad_Config.instance();
                   matRad_cfg.dispError('All objectives have been optimized already in this list');
               end
            

            end

        end

        function cst = getCst(this, list, originalCst, objective, currPriority)
            % This function should build the cst accoriding to the priority

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
            for pIdx=1:currPriority-1

                currObjective = list.getObjectivesWithPriority(pIdx);

                for i=1:numel(currObjective)

                    if ~isempty(currObjective(i).optimValue) % If it's empty we have a problem
                       maxObj = roundObjective(currObjective(i).optimValue*this.slack,8); % Thi is where Phase1 differes from nominal
                    end

                    cstObj            = currObjective(i).objectiveFunction.turnIntoLexicographicConstraint(maxObj);
                    cstObj.quantity   = currObjective(i).quantity;
                    cstObj.robustness = currObjective(i).robustness;

                    cst{currObjective(i).cstIdx, 6} = [cst{currObjective(i).cstIdx, 6}, {cstObj}];
                end

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