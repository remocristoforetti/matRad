classdef matRad_TotalEffectVariance < OmegaObjectives.matRad_OmegaObjective
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
        function obj = matRad_TotalEffectVariance()
        end
        
        function f = computeTotalVarianceObjective(obj,totVariance,d_i)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            f = 1/numel(d_i) * abs(totVariance - d_i'*d_i);


        end
        
        
        function [g,dExpG] = computeTotalVarianceGradient(obj,totVariance,d_i)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            g = (1/numel(d_i)).*sign(totVariance - d_i'*d_i);
            dExpG = - 2 * (1/numel(d_i))*d_i.*sign(totVariance - d_i'*d_i);

        end
    end
end

