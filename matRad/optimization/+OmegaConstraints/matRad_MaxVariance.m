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
        quantity;
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
            s = struct@DoseConstraints.matRad_DoseConstraint(this);
        end

        function jstruct = getDoseConstraintJacobianStructure(~, n)
            jstruct = ones(n, 1);
        end
        
      
 
        
        %% Calculates the Constraint Function value
        function cMeanVariance = computeVarianceConstraintFunction(~,vTot, ~)
            cMeanVariance = vTot;
        end
        
        %% Calculates the Constraint jacobian
        function cVarianceJacob  = computeVarianceConstraintJacobian(~,dOmega, ~)
            if ~isscalar(dOmega)
                % For older code compatibility
                cVarianceJacob = dOmega;
            else
                cVarianceJacob = 1;
            end
        end

        %% Get bounds
        function cl = lowerBounds(~,~)
            cl = 0;
        end
    
        function cu = upperBounds(this,~)
            cu = this.parameters{1};
        end
    
    end
end