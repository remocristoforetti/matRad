classdef matRad_LexicographicObjective < handle

    properties
        objectiveFunction;
        priority;
        goal;
        cstIdx;
        quantity;
        robustness;
    end

    properties (SetAccess = protected)
        optimValue;
    end

    methods

        function this = matRad_LexicographicObjective(priority,objectiveFunction,cstIdx,quantity, goal, varargin)

            p = inputParser();

            addRequired(p, 'priority',     @(x) isnumeric(x));
            addRequired(p, 'objectiveFunction', @(x) isa(x, 'matRad_DoseOptimizationFunction') || isa(x,'OmegaObjectives.matRad_OmegaObjective') || isa(x, 'OmegaConstraints.matRad_VarianceConstraint'));
            addRequired(p, 'cstIdx',       @(x) isnumeric(x));
            addRequired(p, 'quantity',     @(x) ischar(x));
            addRequired(p, 'goal',         @(x) isnumeric(x));
            addOptional(p, 'robustness', 'none', @(x) ischar(x));

            parse(p, priority, objectiveFunction,cstIdx,quantity, goal, varargin{:});

            this.objectiveFunction = p.Results.objectiveFunction;
            this.priority          = p.Results.priority;
            this.cstIdx            = p.Results.cstIdx;
            this.quantity          = p.Results.quantity;
            this.goal              = p.Results.goal;
 
            this.robustness        = p.Results.robustness;

        end

    end


    methods
        function setOptimValue(this, value)

            this.optimValue = value;
        end

        function setGoalValue(this, value)
            this.goal = value;
        end

        function resetOptimValue(this)
            this.optimValue = [];
        end
    end
end