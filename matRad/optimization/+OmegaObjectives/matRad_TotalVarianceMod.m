classdef matRad_TotalVarianceMod < OmegaObjectives.matRad_OmegaObjective
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
        function obj = matRad_TotalVarianceMod()
        end
        
        function f = computeTotalVarianceObjective(obj,totVariance,d_i)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            f = 1/numel(d_i) * abs(totVariance);
        end
        
        
        function [g] = computeTotalVarianceGradient(obj,totalVariance,d_i)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            % g = 1/numel(d_i)*sign(totalVariance);
            g = 1/numel(d_i)*sign(totalVariance);
        end
    end
end

