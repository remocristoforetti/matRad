classdef matRad_TotalVariancePositive < OmegaObjectives.matRad_OmegaObjective
    %MATRAD_MINTOTALVARIANCE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name = 'Min. Abs Total Variance';
        parameterNames = {};
        parameterTypes = {};
    end
    
    properties
        parameters = {0};
        penalty = 1;
        robustness = 'PROB';
    end
    
    methods
        function obj = matRad_TotalVariancePositive()
        end
        
        function f = computeTotalVarianceObjective(obj,totVariance,d_i)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if totVariance>0
                f = 1/numel(d_i) * totVariance;
            else
                f=0;
            end
        end
        
        
        function [g] = computeTotalVarianceGradient(obj,totalVariance,d_i)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            % g = 1/numel(d_i)*sign(totalVariance);
            if totalVariance>0
                g = (1/numel(d_i))*totalVariance;
            else
                g = 0;
            end
        end
    end
end

