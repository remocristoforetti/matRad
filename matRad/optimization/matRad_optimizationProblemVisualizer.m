classdef matRad_optimizationProblemVisualizer < handle
    properties
        hAx_costFunction;
        hAx_costraintFunction;
        hFig;

        hFig_const;

        hPlots_costFunctions;
        hPlots_costraintFunctions;

        active = false;
    end


    properties
        isOpen = false;
        plotFailed = false;
        isOpenConst = false;
        plotConstraintFailed = false;
    end

    properties (SetAccess = protected)
        nObjFunctions;
        cst_structures;
        data;


        leg;
        iter = 0;
        iter_const = 0;
        leg_constr;
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
                this.data.constraintFunctions = [];

                structIdx = find(~cellfun(@isempty, cst(:,6)))';

                for i=structIdx
                    for j=1:numel(cst{i,6})
                        objective = cst{i,6}{j};

                        if isa(objective, 'DoseObjectives.matRad_DoseObjective') || isa(objective, 'OmegaObjectives.matRad_OmegaObjective')
                            
                            objFuncStruct.name = objective.name;
                            objFuncStruct.values = [];
                            
                            this.data.objectiveFunctions = [this.data.objectiveFunctions, objFuncStruct];
                            this.leg = [this.leg, {[cst{i,2}, ' ', objFuncStruct.name]}];
                        elseif isa(objective, 'DoseConstraints.matRad_DoseConstarint') || isa(objective, 'OmegaConstraints.matRad_VarianceConstraint')
                            
                            costFuncStruct.name = objective.name;
                            costFuncStruct.values = [];

                            
                            if isa(objective, 'DoseConstraints.matRad_DoseConstarint')
                                costFuncStruct.ubound = cst{i,6}{j}.upperBounds(numel(cst{i,4}{1}));
                                costFuncStruct.lbound = cst{i,6}{j}.lowerBounds(numel(cst{i,4}{1}));
                            else
                                costFuncStruct.ubound = cst{i,6}{j}.upperBounds;
                                costFuncStruct.lbound = cst{i,6}{j}.lowerBounds;
                            end
                            this.data.constraintFunctions = [this.data.constraintFunctions, costFuncStruct];
                            this.leg_constr= [this.leg_constr, {[cst{i,2}, ' ', costFuncStruct.name]}];

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

        function updateDataConstraints(this, constraintFunctionValues, uBound, lBound)
                this.iter_const = this.iter_const +1;
                if numel(constraintFunctionValues) == numel(this.data.constraintFunctions)
                    for i=1:numel(this.data.constraintFunctions)
                        this.data.constraintFunctions(i).values(this.iter_const) = constraintFunctionValues(i);
                    
                        if isempty(this.data.constraintFunctions(i).ubound)
                            this.data.constraintFunctions(i).ubound = uBound;
                        end
    
                        if isempty(this.data.constraintFunctions(i).lbound)
                            this.data.constraintFunctions(i).lbound = lBound;
                        end
                    end
                   
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

        function updatePlotConstraints(this)
            
            matRad_cfg = MatRad_Config.instance();
            
            if ~this.plotConstraintFailed
                try
                    this.plotFunctionConstr();
                catch
                    matRad_cfg.dispWarning('Single Objective Function plotting failed and thus disabled.');
                    this.plotConstraintFailed = true;
                    this.isOpenConst = false;
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
                    %set(curr_hAx,'YScale','log');
                    title(curr_hAx,'Progress of Optimization','LineWidth',defaultFontSize);
                    xlabel(curr_hAx,'# iterations','Fontsize',defaultFontSize),ylabel(curr_hAx,'objective function value','Fontsize',defaultFontSize);
                    
                    %Create plot handle and link to data for faster update
                    for i=1:numel(this.data.objectiveFunctions)
                        y = this.data.objectiveFunctions(i).values;
                        hPlot = plot(curr_hAx,x,y,'x--','LineWidth',0.5,'XDataSource','x','YDataSource','y');
                        this.hPlots_costFunctions = [this.hPlots_costFunctions, hPlot];
                    end
                    y = this.data.totFValues;
                    hPlot = plot(curr_hAx,x,y,'x--','LineWidth',0.5,'XDataSource','x','YDataSource','y');
                    this.hPlots_costFunctions = [this.hPlots_costFunctions, hPlot];
    
                    legend(curr_hAx, [this.leg, 'totalFunction']);
                    
                    this.hAx_costFunction = curr_hAx;
    
                    this.hFig = curr_hFig;
                    this.isOpen = true;
                else
    
                    curr_hFig = get(this.hAx_costFunction,'Parent');
                    curr_hAx = this.hAx_costFunction;
                    for plotIdx = 1:numel(this.data.objectiveFunctions)
            
                        y = this.data.objectiveFunctions(plotIdx).values;
                        
                        curr_hPlot = this.hPlots_costFunctions(plotIdx);
                        refreshdata(curr_hPlot, 'caller');
        
                    end
                    y = this.data.totFValues;
                    curr_hPlot = this.hPlots_costFunctions(end);
                    refreshdata(curr_hPlot, 'caller');

                end
    
                drawnow;
                figure(curr_hFig);
                %movegui(curr_hFig, 'southwest');
            
        end


        function plotFunctionConstr(this)

            
                x = [1:this.iter_const];
                if ~this.isOpenConst
                
                    curr_hFig = figure('Name','Progress of single constraints','NumberTitle','off');

                    curr_hAx  = axes(curr_hFig);
        
                    hold(curr_hAx,'on');
                    grid(curr_hAx,'on');
                    grid(curr_hAx,'minor');
                       
                    %Set up the axes scaling & labels
                    defaultFontSize = 14;
                    %title(curr_hAx,'','LineWidth',defaultFontSize);
                    xlabel(curr_hAx,'# iterations','Fontsize',defaultFontSize),ylabel(curr_hAx,'constraint function value','Fontsize',defaultFontSize);
                    
                    %Create plot handle and link to data for faster update
                    for i=1:numel(this.data.constraintFunctions)
                        if ~isempty(this.data.constraintFunctions(i).values)
                            y = this.data.constraintFunctions(i).values;
                            hPlot = plot(curr_hAx,x,y,'x--','LineWidth',0.5,'XDataSource','x','YDataSource','y');
                            
                            yline(this.data.constraintFunctions(i).ubound);
                            yline(this.data.constraintFunctions(i).lbound);
    
                            % y = this.data.constraintFunctions(i).ubound;
                            % plot(curr_hAx,x,y,'.--','LineWidth',0.5,'XDataSource','x','YDataSource','y');
                            % y = this.data.constraintFunctions(i).lbound;
                            % plot(curr_hAx,x,y,'.--','LineWidth',0.5,'XDataSource','x','YDataSource','y');
    
                            this.hPlots_costraintFunctions = [this.hPlots_costraintFunctions, hPlot];
                        end
                    end
                    
                    this.hAx_costraintFunction = curr_hAx;
    
                    this.hFig_const = curr_hFig;
                    this.isOpenConst = true;
                else
                    curr_hFig = get(this.hAx_costraintFunction,'Parent');    
                    curr_hAx = this.hAx_costraintFunction;
                    for plotIdx = 1:numel(this.data.constraintFunctions)
            
                        y = this.data.constraintFunctions(plotIdx).values;
                        
                        curr_hPlot = this.hPlots_costraintFunctions(plotIdx);
                        refreshdata(curr_hPlot, 'caller');
        
                    end
    
                end
    
                drawnow;
                figure(curr_hFig);
                %movegui(curr_hFig, 'southwest');
            end


        %% Get
        function isOpen = get.isOpen(this)
            if ~isempty(this.hFig)
                isOpen = true;
            else
                isOpen = false;
            end
        end

        function isOpen = get.isOpenConst(this)
            if ~isempty(this.hFig_const)
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

                    this.isOpenConst = false;
                    this.plotConstraintFailed = true;
                end
            end
        end
    end
end