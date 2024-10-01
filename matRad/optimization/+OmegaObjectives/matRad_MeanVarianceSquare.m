classdef matRad_MeanVarianceSquare < OmegaObjectives.matRad_OmegaObjective
    %MATRAD_MINTOTALVARIANCE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name = 'Min. Total Variance';
        parameterNames = {};
        parameterTypes = {};
    end
    
    properties
        parameters = {0};
        penalty = 1;
        robustness = 'PROB';
    end
    
    methods
        function obj = matRad_MeanVarianceSquare()
        end
        
        function f = computeTotalVarianceObjective(obj,totVariance,d_i)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            f = totVariance^2;
        end
        
        
        function g = computeTotalVarianceGradient(obj,totVariance,d_i)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            g = 2 * totVariance;
        end
    end
end

