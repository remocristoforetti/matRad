classdef matRad_EUD_SqOver < DoseObjectives.matRad_DoseObjective
% matRad_EUD Implements a penalized equivalent uniform dose objective
%   See matRad_DoseObjective for interface description
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties (Constant)
        name = 'EUD';
        parameterNames = {'EUD^{ref}', 'k', 'd^{max}'};
        parameterTypes = {'dose','numeric', 'dose'};
    end
    
    properties
        parameters = {0, 3.5, 0};
        penalty;
    end
    
    methods
        function obj = matRad_EUD_SqOver(penalty,eudRef, eudExponent, doseRef)
            %If we have a struct in first argument
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            obj@DoseObjectives.matRad_DoseObjective(inputStruct);
            
            %now handle initialization from other parameters
            if ~initFromStruct
                if exist('penalty', 'var') && ~isempty(penalty)
                    obj.penalty = penalty;
                end

                if exist('eudRef', 'var') && ~isempty(eudRef)
                    obj.parameters{1} = eudRef;
                end

                if exist('eudExponent', 'var') && ~isempty(eudExponent)
                    obj.parameters{2} = eudExponent;
                end
                
                if exist('doseRef', 'var') && ~isempty(doseRef)
                    obj.parameters{3} = doseRef;
                end
            end
        end
        
        %% Calculates the Objective Function value
        function fDose = computeDoseObjectiveFunction(obj,dose)
            % get exponent for EUD
            k = obj.parameters{2};
            
            % calculate power sum

            if any(dose<0)
                matRad_cfg = MatRad_Config.instance();
                stringWarning = sprintf('Negative dose values detected:%s, setting to zero', dose(dose<0));
                matRad_cfg.dispWarning(stringWarning);
                %dose = dose(dose>0);
                dose(dose<0) = 0;
            end
            powersum = sum(dose.^k);
            
            %Calculate objective
            
            %This check is not needed since dose is always positive
            %if powersum > 0
            fDose = (nthroot(powersum/numel(dose),k) - obj.parameters{1})^2;


            overdose = dose - obj.parameters{3};
            
            % apply positive operator
            overdose(overdose<0) = 0;
            
            % claculate objective function
            fDose = fDose + 1/numel(dose) * (overdose'*overdose);
            %end
        end
        
        %% Calculates the Objective Function gradient
        function fDoseGrad  = computeDoseObjectiveGradient(obj,dose)
            
            overdose = dose - obj.parameters{3};
            
            % apply positive operator
            overdose(overdose<0) = 0;

            % get exponent for EUD
            k = obj.parameters{2};
            
            %numerical stability
            dose(dose == 0) = 0.001;
            
            if any(dose<0)
                matRad_cfg = MatRad_Config.instance();
                stringWarning = sprintf('Negative dose values detected:%s, setting to zero', dose(dose<0));
                matRad_cfg.dispWarning(stringWarning);
                %dose = dose(dose>0);
                dose(dose<0) = 0.001;
            end
            
            % calculate power sum
            powersum = sum(dose.^k);
                        
            
            %This check is not needed since dose is always positive
            %if powersum > 0
            
            %derivatives = nthroot(1/numel(dose),k) * powersum^((1-k)/k) * (dose.^(k-1));
            fDoseGrad = 2 * nthroot(1/numel(dose),k) * powersum^((1-k)/k) * (dose.^(k-1)) .* (nthroot(powersum/numel(dose),k) - obj.parameters{1});
            %end

            

            
            % calculate delta
            fDoseGrad = fDoseGrad + 2 * 1/numel(dose) * overdose;
            if any(~isfinite(fDoseGrad)) % check for inf and nan for numerical stability
                error(['EUD computation failed. Reduce exponent to resolve numerical problems.']);
            end
        end

        function constr = turnIntoLexicographicConstraint(obj,goal)
            objective = DoseObjectives.matRad_EUD_SqOver(100,obj.parameters{1,1},obj.parameters{1,2},obj.parameters{1,3});
            objective.quantity = obj.quantity;
            objective.robustness = obj.robustness;
            constr = DoseConstraints.matRad_DoseConstraintFromObjective(objective,goal);
            
            constr.quantity = obj.quantity;
            constr.robustness = obj.robustness;

        end

    end
    
    methods (Static)
        function newGoalValue = adaptGoalToFraction(goalValue,numOfFractions)
            newGoalValue = goalValue/numOfFractions;
        end
        
    end
    
end

