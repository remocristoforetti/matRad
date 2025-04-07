classdef matRad_VisualizationManager < handle

    properties
        hFig;
        modules;
        nModules;
    
        dataQuery;
    end

    properties (SetAccess=protected)
        isOpen = false;
        plotFailed = false;
    end

    methods
        
        function this = matRad_VisualizationManager()
            
            % valid = arrayfun(@(x) isa(x, 'matRad_VisualizationModule'));
            % 
            % if all(valid)
            %     this.modules  = modules;
            %     this.nModules = numel(modules);
            % end
            
            this.hFig = figure('Name','Visualization Manager','NumberTitle','off','Visible','off');
            tiledlayout('flow');
            this.hFig.Position(3:4) = [900,900];
            this.hFig.Position(1:2) = this.hFig.Position(1:2)- 500;
        end

        function updatePlot(this)
            
            matRad_cfg = MatRad_Config.instance();
            
            if ~this.plotFailed
                
                try
                    this.hFig.Visible = 'on';
                    for curModule=1:numel(this.modules)
                        this.modules{curModule}.plotData();
                    end
                catch
                    matRad_cfg.dispWarning('Visualization manager disabled.');
                    this.plotFailed = true;
                    this.isOpen = false;
                    close(this.hFig);
                end

            else
                % if isvalid(this.hFig)
                %     this.hFig
                % end
            end

        end

        function addModule(this, module)
            if isa(module, 'matRad_VisualizationModule')
                this.modules = [this.modules, {module}];
                this.modules{end}.hAx = nexttile(this.hFig.Children);
            end
            % Else...
        end

        function removeModule(this, module)

        end

        function delete(this)
            delete(this.hFig);
        end
    end
end 