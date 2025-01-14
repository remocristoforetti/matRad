classdef matRad_TotalVariance < OmegaObjectives.matRad_OmegaObjective
    %MATRAD_MINTOTALVARIANCE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name = 'Min. Total Variance';
        parameterNames = {};
        parameterTypes = {};
    end
    
    properties
        parameters = {};
        penalty;
        robustness = 'PROB';
    end
    
    methods
        function obj = matRad_TotalVariance(penalty)

            if exist('penalty', 'var') && ~isempty(penalty)
                obj.penalty = penalty;
            else
                obj.penalty = 1;
            end
        end
        
        function f = computeTotalVarianceObjective(obj,totVariance,nVoxels)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            f = (10)*1/nVoxels * totVariance;
        end
        
        
        function g = computeTotalVarianceGradient(obj,~,nVoxels)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            g = (10)*1/nVoxels;
        end

        function constr = turnIntoLexicographicConstraint(obj,goal)
            if goal < 5e-4
                goal = 5e-4*1.03;
            end
            
            objective = OmegaObjectives.matRad_TotalVariance();
            objective.quantity = obj.quantity;
            objective.robustness = obj.robustness;
            
            constr = OmegaConstraints.matRad_OmegaConstraintFromObjective(objective,goal);
            constr.quantity = obj.quantity;
            constr.robustness = obj.robustness;
        end
    end
end

