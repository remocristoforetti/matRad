classdef matRad_MeanVariance < OmegaObjectives.matRad_OmegaObjective
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
        function obj = matRad_MeanVariance()
        end
        
        function f = computeTotalVarianceObjective(obj,totVariance,d_i)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if totVariance>0
                f = totVariance;
            else
                f = 0;
            end
        end
        
        
        function g = computeTotalVarianceGradient(obj,totVariance,d_i)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if totVariance>0
                g = 1;
            else
                g = 0;
            end
        end
    end
end

