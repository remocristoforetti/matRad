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
        penalty = 1;
        robustness = 'PROB';
    end
    
    methods
        function obj = matRad_TotalVariance()
        end
        
        function f = computeTotalVarianceObjective(obj,totVariance,d_i)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            f = 1/numel(d_i) * totVariance;
        end
        
        
        function [g,dExpG] = computeTotalVarianceGradient(obj,~,d_i)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            g = 1/numel(d_i);
            dExpG = 0;
        end
    end
end

