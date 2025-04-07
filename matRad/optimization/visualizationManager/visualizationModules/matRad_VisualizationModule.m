classdef matRad_VisualizationModule < handle

    properties
        hAx;
    end

    methods
        function this = matRad_VisualizationModule()

        end

        function updateData(this)
            % To be implemented by subclasses
        end

        function plotData(this)

        end
    end
end