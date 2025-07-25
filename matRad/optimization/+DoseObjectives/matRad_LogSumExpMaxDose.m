classdef matRad_LogSumExpMaxDose < DoseObjectives.matRad_DoseObjective
% matRad_MeanDose Implements a penalized MeanDose objective
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
        name = 'Max Dose';
        parameterNames = {'d^{ref}','f_{diff}'}; %When optimizing to a reference, one might consider using a quadratic relationship with a non-linear optimizer
        parameterTypes = {'dose',{'Linear','Quadratic'}};
    end
    
    properties
        parameters = {0,1};
        penalty = 1;
        epsilon = 1e-3;
        % epsilon = 1e-4;
    end
    
    methods 
        function obj = matRad_LogSumExpMaxDose(penalty,dMaxRef,fDiff)
           
            % if we have a struct in first argument
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            obj@DoseObjectives.matRad_DoseObjective(inputStruct);
            
            if ~initFromStruct
                if nargin < 3 || ~ischar(fDiff)
                    fDiff = 'Linear';
                end
                
                fDiffIx = find(strcmp(fDiff,obj.parameterTypes{2}));
                
                if isempty(fDiffIx) || numel(fDiffIx) > 1
                    fDiffIx = 1;
                    matRad_cfg = MatRad_Config.instance();                    
                    matRad_cfg.dispWarning('Mean dose difference function can only be %s! Using %s difference.', strjoin(obj.parameterTypes{2},' or '), obj.parameterTypes{2}{fDiffIx});
                end
                
                obj.parameters{2} = fDiffIx;


                if nargin >= 2 && isscalar(dMaxRef)
                    obj.parameters{1} = dMaxRef;
                end

                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end

            %% Downwards compatability / set default values
            %TODO: maybe move into set method for parameters
            if numel(obj.parameters) < 1
                obj.parameters{1} = 0;
            end

            if numel(obj.parameters) < 2
                obj.parameters{2} = 1;
            end
            
        end       
        
        %% Calculates the Objective Function value
        function fDose = computeDoseObjectiveFunction(obj,dose)
            switch obj.parameters{2}
                case 1
                    fDose =  obj.objectiveLinearDiff(dose);
                case 2
                    fDose =  obj.objectiveQuadraticDiff(dose);
                otherwise
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Invalid setting for %s in Mean Dose Objective!',obj.parameterNames{2});  
            end
        end
        
        %% Calculates the Objective Function gradient
        function fDoseGrad   = computeDoseObjectiveGradient(obj,dose)
            switch obj.parameters{2}
                case 1
                    fDoseGrad = obj.gradientLinearDiff(dose);
                case 2
                    fDoseGrad = obj.gradientQuadraticDiff(dose);
                otherwise
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Invalid setting for %s in Mean Dose Objective!',obj.parameterNames{2});  
            end
        end

        function constr = turnIntoLexicographicConstraint(obj,goal)
            objective = DoseObjectives.matRad_MaxDose(100,obj.parameters{1},obj.parameters{2});
            objective.quantity = obj.quantity;
            objective.robustness = obj.robustness;
            constr = DoseConstraints.matRad_DoseConstraintFromObjective(objective,goal);
            constr.quantity = obj.quantity;

            constr.robustness = obj.robustness;

        end

    end

    methods (Access = protected)
        function fDose = objectiveQuadraticDiff(obj,dose)
            fDose = (max(dose(:)) - obj.parameters{1})^2;
        end

        function fDoseGrad = gradientQuadraticDiff(obj,dose)
            fDoseGrad = 2*(max(dose(:))-obj.parameters{1}) * ones(size(dose(:)))/numel(dose);
        end

        function fDose = objectiveLinearDiff(obj,dose)
            dose_max = max(dose);
            modEpsilon = (obj.epsilon)*dose_max;
            fDose = dose_max + modEpsilon * log( sum(exp((dose - dose_max)/modEpsilon)));
        end

        function fDoseGrad = gradientLinearDiff(obj,dose)
            [max_dose,maxDoseIdx] = max(dose);
            modEpsilon = (obj.epsilon)*max_dose;

            fDoseGrad(:,1) = exp( (dose-max_dose)/modEpsilon );
            fDoseGrad(:,1) = fDoseGrad(:,1)/sum(fDoseGrad(:,1));
            % vExp = exp( (dose-max_dose)/modEpsilon);
            % vSumExp = sum(vExp);
            % 
            % fDoseGrad(:,1) = (modEpsilon/vSumExp) * vExp;
            % fDoseGrad(maxDoseIdx,1) = 1 + ...
            %                            (modEpsilon/max_dose)*log(vSumExp) + ...
            %                            (1/vSumExp)*(1 - vSumExp - (modEpsilon/max_dose)*sum(vExp.*((dose-max_dose)/modEpsilon)));

        end
    end

    methods (Static)
        function newGoalValue = adaptGoalToFraction(goalValue,numOfFractions)
            newGoalValue = goalValue/numOfFractions;
        end
    end
end

