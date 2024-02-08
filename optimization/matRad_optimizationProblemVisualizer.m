classdef matRad_optimizationProblemVisualizer < handle
    properties
        hAx;
        hFig;
        hPlots;

        active = false;
    end


    properties
        isOpen = false;
        plotFailed = false;
    end

    properties (SetAccess = protected)
        nObjFunctions;
        cst_structures;
        data;


        leg;
        iter = 0;
    end


    methods
        function this = matRad_optimizationProblemVisualizer(cst)

            matRad_cfg = MatRad_Config.instance();

            % if isnumeric(nObjFunctions) && ~isempty(nObjFunctions)
            %     this.nObjFunctions = nObjFunctions;
            % else
            %     matRad_cfg.dispWarning('Visualizer: Empty or non valid number of ObjectiveFunctions, setting to default')
            %     this.nObjFunctions = 1;
            % end

            if ~isempty(cst)
                
                this.data.objectiveFunctions = [];
                structIdx = find(~cellfun(@isempty, cst(:,6)))';

                for i=structIdx
                    for j=1:numel(cst{i,6})
                        objective = cst{i,6}{j};

                        if isa(objective, 'DoseObjectives.matRad_DoseObjective') || isa(objective, 'OmegaObjectives.matRad_OmegaObjective')
                            
                            objFuncStruct.name = objective.name;
                            objFuncStruct.values = [];
                            
                            this.data.objectiveFunctions = [this.data.objectiveFunctions, objFuncStruct];
                            this.leg = [this.leg, {[cst{i,2}, ' ', objFuncStruct.name]}];
                        end
                    end
                end
            end
        end

        function updateData(this, objFunctionValues, totF)

            this.iter = this.iter +1;

            if numel(objFunctionValues) == numel(this.data.objectiveFunctions)
                for i=1:numel(this.data.objectiveFunctions)
                    this.data.objectiveFunctions(i).values(this.iter) = objFunctionValues(i);
                
                end
                this.data.totFValues(this.iter) = totF;
            end
        end

        
        function updatePlot(this)
            
            matRad_cfg = MatRad_Config.instance();
            if ~this.plotFailed
                try 

                    this.plotFunction();
                catch
                    matRad_cfg.dispWarning('Single Objective Function plotting failed and thus disabled.');
                    this.plotFailed = true;
                    this.isOpen = false;
                end

            end

        end

        function plotFunction(this)

            x = [1:this.iter];

            if ~this.isOpen
                curr_hFig = figure('Name','Progress of single objectives','NumberTitle','off');
                curr_hAx  = axes(curr_hFig);

                hold(curr_hAx,'on');
                grid(curr_hAx,'on');
                grid(curr_hAx,'minor');
                set(curr_hAx,'YScale','log');

                %Set up the axes scaling & labels
                defaultFontSize = 14;
                set(curr_hAx,'YScale','log');
                title(curr_hAx,'Progress of Optimization','LineWidth',defaultFontSize);
                xlabel(curr_hAx,'# iterations','Fontsize',defaultFontSize),ylabel(curr_hAx,'objective function value','Fontsize',defaultFontSize);
                
                %Create plot handle and link to data for faster update
                for i=1:numel(this.data.objectiveFunctions)
                    y = this.data.objectiveFunctions(i).values;
                    hPlot = plot(curr_hAx,x,y,'x--','LineWidth',0.5,'XDataSource','x','YDataSource','y');
                    this.hPlots = [this.hPlots, hPlot];
                end
                y = this.data.totFValues;
                hPlot = plot(curr_hAx,x,y,'x--','LineWidth',0.5,'XDataSource','x','YDataSource','y');
                this.hPlots = [this.hPlots, hPlot];

                legend(curr_hAx, [this.leg, 'totalFunction']);
                this.hFig = curr_hFig;
                this.hAx = curr_hAx;
                this.isOpen = true;
            else

                curr_hFig = get(this.hAx,'Parent');
                curr_hAx = this.hAx;
                for plotIdx = 1:numel(this.data.objectiveFunctions)
        
                    y = this.data.objectiveFunctions(plotIdx).values;
                    
                    curr_hPlot = this.hPlots(plotIdx);
                    refreshdata(curr_hPlot, 'caller');
    
                end
                y = this.data.totFValues;

                curr_hPlot = this.hPlots(end);
                refreshdata(curr_hPlot, 'caller');
            end

            drawnow;
            figure(curr_hFig);
            movegui(curr_hFig, 'southwest');
        end

        %% Get
        function isOpen = get.isOpen(this)
            if ~isempty(this.hFig)
                isOpen = true;
            else
                isOpen = false;
            end
        end

        %% Set 
        function set.active(this, boolValue)
            valid = islogical(boolValue);

            if valid
                this.active = boolValue;
                if ~boolValue
                    this.isOpen = false;
                    this.plotFailed = true;
                end
            end
        end
    end
end