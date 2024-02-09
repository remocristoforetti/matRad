classdef (Abstract) matRad_Robustness < handle
    properties
        useScen; % scenario idx on the dij.physicalDose matrix
    end


    methods
        
        function this = matRad_Robustness()
            % if nargin>0 && ~isempty()
            %     this.assignScenarios();
            % end
        end

        function assignScenarios(this,multScen)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispWarning('Robustness subclass needs to implement assignScenario method');

        end


    end

end