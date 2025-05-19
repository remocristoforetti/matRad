classdef matRad_ObjectiveFunctionVisModule < matRad_VisualizationModule

    properties
        data;
        leg;
        leg_constr;

        dataRequest = {'singleObjective', 'f'};
    end

    properties (Hidden)
        defaultFontSize = 14;
    end

    methods
        function this = matRad_ObjectiveFunctionVisModule(cst)
            matRad_cfg = MatRad_Config.instance();

            if ~isempty(cst)
                
                this.data.objectiveFunctions = [];
                this.data.constraintFunctions = [];

                structIdx = find(~cellfun(@isempty, cst(:,6)))';

                for i=structIdx
                    for j=1:numel(cst{i,6})
                        objective = cst{i,6}{j};

                        if isstruct(objective)
                            objective = matRad_DoseOptimizationFunction.createInstanceFromStruct(objective);
                        end

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
            this.data.iter = 0;
        end

        function updateData(this, dataInput)

            this.updateDataObjectives(dataInput.singleObjective, dataInput.f);

        end
        function updateDataObjectives(this, objFunctionValues, totF)

            this.data.iter = this.data.iter +1;
            if numel(objFunctionValues) == numel(this.data.objectiveFunctions)
                for i=1:numel(this.data.objectiveFunctions)
                    this.data.objectiveFunctions(i).values(this.data.iter) = objFunctionValues(i);
                
                end
                if exist('totF', 'var')
                    this.data.totFValues(this.data.iter) = totF;
                else
                    this.data.totFValues(this.data.iter) = sum([this.data.objectiveFunctions(:).values(this.data.iter)]);
                end
            end
        end

        function updateDataConstraints(this, constraintFunctionValues, uBound, lBound)
            
                this.data.iter_const = this.data.iter_const +1;
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

        function plotData(this)

                x = [1:this.data.iter];

                title('Single objectives');
    
                hold(this.hAx,'on');
                grid(this.hAx,'on');
                grid(this.hAx,'minor');
                set(this.hAx,'YScale','log');
    
                %Set up the axes scaling & labels
                this.defaultFontSize = 14;

                title(this.hAx,'Progress of Optimization','LineWidth',this.defaultFontSize);
                xlabel(this.hAx,'# iterations','Fontsize',this.defaultFontSize),ylabel(this.hAx,'objective function value','Fontsize',this.defaultFontSize);
                    
                c = colororder;
                %Create plot handle and link to data for faster update
                for i=1:numel(this.data.objectiveFunctions)
                    y = this.data.objectiveFunctions(i).values;
                    hPlot = plot(this.hAx,x,y,'x--','LineWidth',0.5,'XDataSource','x','YDataSource','y', 'Color', c(i,:));
                    %this.hPlots_costFunctions = [this.hPlots_costFunctions, hPlot];
                end
                y = this.data.totFValues;
                hPlot = plot(this.hAx,x,y,'x--','LineWidth',0.5,'XDataSource','x','YDataSource','y', 'Color', c(i+1,:));
                %this.hPlots_costFunctions = [this.hPlots_costFunctions, hPlot];

                legend(this.hAx, [this.leg, 'totalFunction']);
                    

                    %%% Dose distribution
                    % subplot(1,2,2);
                    % curr_hAx_dist = gca();
                    % 
                    % doseGrid = this.distibutionProperties.doseGrid;
                    % ctGrid   = this.distibutionProperties.ct;
                    % 
                    % dDoseGrid = reshape(this.d.(this.distibutionProperties.quantity){1}, doseGrid.dimensions);
                    % dCtGrid = matRad_interp3(doseGrid.x,doseGrid.y',doseGrid.z, ...
                    %                          dDoseGrid, ...
                    %                          ctGrid.x,ctGrid.y',ctGrid.z,'linear',0);
                    % 
                    % matRad_plotSliceWrapper(curr_hAx_dist, this.distibutionProperties.ct, this.distibutionProperties.cst, 1, dCtGrid, this.distibutionProperties.plane, this.distibutionProperties.slice, [],[],[],[],[],[],ones(1,size(this.distibutionProperties.cst,1)));
                    % 
                    %this.hAx_costFunction = curr_hAx;
                    %this.hAx_distribution = curr_hAx_dist;
                    %this.hFig = curr_hFig;
                    %this.isOpen = true;
                drawnow;
                %figure(curr_hFig);
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
    end


end