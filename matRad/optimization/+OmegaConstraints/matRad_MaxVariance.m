classdef matRad_MaxVariance < OmegaConstraints.matRad_VarianceConstraint
    % matRad_MinMaxDose Implements a MinMaxDose constraint
    %   See matRad_DoseConstraint for interface description
    %
    % use log sum exp approximation, see appendix A in
    % http://scitation.aip.org/content/aapm/journal/medphys/41/8/10.1118/1.4883837
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
        name = 'Max Varaince constraint';
        parameterNames = {'var^{max}'};
        parameterTypes = {'maxVariance'};
    end
    
    properties
        parameters = {30};
        robustness;
    end
    
    methods
        function this = matRad_MaxVariance(maxVaraince)
            if exist('maxVaraince', 'var') && ~isempty(maxVaraince)
                this.parameters = {maxVaraince};
            end
        end
        
        %Overloads the struct function to add constraint specific
        %parameters
        function s = struct(this)
 %           s = struct@DoseConstraints.matRad_DoseConstraint(this);
             s = struct@OmegaConstraints.matRad_VarianceConstraint(this);
        end

        function jstruct = getDoseConstraintJacobianStructure(~, n)
            jstruct = ones(n, 1);
        end
        
      
 
        
        %% Calculates the Constraint Function value
        function cMeanVariance = computeVarianceConstraintFunction(~,vTot, nVoxels)
            cMeanVariance = vTot/nVoxels;
        end
        
        %% Calculates the Constraint jacobian
        function cVarainceJacob  = computeVarianceConstraintJacobian(~,dOmega, nVoxels)
            cVarainceJacob = (2/nVoxels) * dOmega;
        end

        %% Get bounds
        function cl = lowerBounds(~)
            cl = 0;
        end
    
        function cu = upperBounds(this)
            cu = this.parameters{1};
        end
    
    end
    
    methods (Access = private)
        % LogSumExp Approximation
        function cDose = computeDoseConstraintFunctionLogSumExp(this,dose)
            dose_min = min(dose);
            dose_max = max(dose);
            
            %Validate parameters
            if this.parameters{1} <= 0 && isinf(this.parameters{2}) %Constraint doesn't make sense (min = 0 & max = Inf)
                cDose = [];
            elseif this.parameters{2} == Inf %Only min dose
                cDose = dose_min - this.epsilon * log( sum(exp((dose_min - dose)/this.epsilon)));
            elseif this.parameters{1} <= 0 %Only max dose
                cDose = dose_max + this.epsilon * log( sum(exp((dose - dose_max)/this.epsilon)));
            else %both are set sensible
                cDose(2,1) = dose_max + this.epsilon * log( sum(exp((dose - dose_max)/this.epsilon)));
                cDose(1,1) = dose_min - this.epsilon * log( sum(exp((dose_min - dose)/this.epsilon)));
            end
            
        end
        function cDoseJacob  = computeDoseConstraintJacobianLogSumExp(this,dose)
            %Validate parameters
            if this.parameters{1} <= 0 && isinf(this.parameters{2}) %Constraint doesn't make sense (min = 0 & max = Inf)
                cDoseJacob = [];
            elseif this.parameters{2} == Inf %Only min dose
                cDoseJacob(:,1) = exp( (min(dose)-dose)/this.epsilon );
                cDoseJacob(:,1) = cDoseJacob(:,1)/sum(cDoseJacob(:,1));
            elseif this.parameters{1} <= 0 %Only max dose
                cDoseJacob(:,1) = exp( (dose-max(dose))/this.epsilon );
                cDoseJacob(:,1) = cDoseJacob(:,1)/sum(cDoseJacob(:,1));
            else %both are set sensible
                cDoseJacob(:,1) = exp( (min(dose)-dose)/this.epsilon );
                cDoseJacob(:,1) = cDoseJacob(:,1)/sum(cDoseJacob(:,1));
                
                cDoseJacob(:,2) = exp( (dose-max(dose))/this.epsilon );
                cDoseJacob(:,2) = cDoseJacob(:,2)/sum(cDoseJacob(:,2));
            end
            
            
        end
        
        %Exact voxel-wise
        function cDose = computeDoseConstraintFunctionVoxelwise(this,dose)
            cDose = dose;
        end
        function cDoseJacob  = computeDoseConstraintJacobianVoxelwise(this,dose)
            cDoseJacob = speye(numel(dose),numel(dose));
        end
    end
    
end


