classdef (Abstract) matRad_OmegaObjective
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Abstract, Access = public)
        penalty %Optimization penalty
        robustness %robustness setting
    end    
    
    methods (Static)
        function rob = availableRobustness()
            rob = {'PROB'}; %By default, no robustness is available
        end 
    end
    
    %These should be abstract methods, however Octave can't parse them. As soon
    %as Octave is able to do this, they should be made abstract again
    methods %(Abstract)
        function f = computeTotalVarianceObjective(obj,totVariance,nVoxels)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            error('Function needs to be implemented!');
        end
        
        
        function g = computeTotalVarianceGradient(obj,totVariance,nVoxels)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            error('Function needs to be implemented!');
        end
    end
    
    methods (Access = public)
       
        % constructor of matRad_DoseObjective
        function obj = matRad_OmegaObjective(varargin)
            %default initialization from struct (parameters & penalty)
            %obj@matRad_DoseOptimizationFunction(varargin{:});
        end
        
        function d = getDoseParameters(obj)
           d = []; 
        end
    end
end

