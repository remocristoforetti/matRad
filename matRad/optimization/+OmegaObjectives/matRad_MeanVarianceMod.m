classdef matRad_MeanVarianceMod < OmegaObjectives.matRad_OmegaObjective
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
        function obj = matRad_MeanVarianceMod()
        end
        
        function f = computeTotalVarianceObjective(obj,totVariance,d_i)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            f = abs(totVariance);
        end
        
        
        function g = computeTotalVarianceGradient(obj,totVariance,d_i)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            g = sign(totVariance);
        end
    end
end

