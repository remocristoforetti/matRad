classdef matRad_ConstraintFunctionVisModule < matRad_VisualizationModule

    properties
        data;
        leg;
        leg_constr;

        dataRequest = {'c'};
    end

    properties (Hidden)
        defaultFontSize = 14;
    end

    methods
        function this = matRad_ConstraintFunctionVisModule(cst)
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

                        if isa(objective, 'DoseConstraints.matRad_DoseConstraint') || isa(objective, 'OmegaConstraints.matRad_VarianceConstraint')
                            
                            if isa(objective, 'DoseConstraints.matRad_DoseConstraintFromObjective')
                                costFuncStruct.name = objective.objective.name;
                                costFuncStruct.values = [];                                
  
                            else
                                costFuncStruct.name = objective.name;
                                costFuncStruct.values = [];
                            end
                            
                            if isa(objective, 'DoseConstraints.matRad_DoseConstraint')
                                costFuncStruct.ubound = min(cst{i,6}{j}.upperBounds(numel(cst{i,4}{1})));
                                costFuncStruct.lbound = max(cst{i,6}{j}.lowerBounds(numel(cst{i,4}{1})));
                            else
                                costFuncStruct.ubound = cst{i,6}{j}.upperBounds;
                                costFuncStruct.lbound = cst{i,6}{j}.lowerBounds;
                            end

                            costFuncStruct.nConstraints = numel(cst{i,6}{j}.upperBounds);
                            this.data.constraintFunctions = [this.data.constraintFunctions, costFuncStruct];
                            this.leg_constr= [this.leg_constr, {[cst{i,2}, ' ', costFuncStruct.name]}];

                            costFuncStruct = [];
                        end
                    end
                end
            end
            
            this.data.iter = 0;
        end

        function updateData(this, dataInput)
            
            this.updateDataObjectives(dataInput);
        end
        
        function updateDataObjectives(this, dataInput)

            this.data.iter = this.data.iter +1;

            if isfield(dataInput, 'c')
                c = dataInput.c;
            else
                c = [];
            end

            if ~isempty(c)
                cCounter = 1;
                for i=1:numel(this.data.constraintFunctions)
                    this.data.constraintFunctions(i).values(this.data.iter,:) = c(cCounter:cCounter+this.data.constraintFunctions(i).nConstraints-1,:);
                    cCounter = cCounter+1;
                end
            end
        end

        % function updateDataConstraints(this, constraintFunctionValues, uBound, lBound)
        % 
        %         this.data.iter_const = this.data.iter_const +1;
        %         if numel(constraintFunctionValues) == numel(this.data.constraintFunctions)
        %             for i=1:numel(this.data.constraintFunctions)
        %                 this.data.constraintFunctions(i).values(this.iter_const) = constraintFunctionValues(i);
        % 
        %                 if isempty(this.data.constraintFunctions(i).ubound)
        %                     this.data.constraintFunctions(i).ubound = uBound;
        %                 end
        % 
        %                 if isempty(this.data.constraintFunctions(i).lbound)
        %                     this.data.constraintFunctions(i).lbound = lBound;
        %                 end
        %             end                 
        %         end
        % 
        % end

        function plotData(this)

                x = [1:this.data.iter];

                title('Single Constraints');
    
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
                for i=1:numel(this.data.constraintFunctions)
                    
                    y = this.data.constraintFunctions(i).values;
                    if ~isempty(y)
                        % hPlot = plot(this.hAx,x,y,'x--','LineWidth',0.5,'XDataSource','x','YDataSource','y', 'Color', c(numel(this.data.objectiveFunctions) + i,:), 'DisplayName', this.data.constraintFunctions(i).name);
                        hPlot = plot(this.hAx,x,y,'x--','LineWidth',0.5,'XDataSource','x','YDataSource','y', 'Color', c(i,:));
                    
                        y1 = yline(this.hAx,this.data.constraintFunctions(i).ubound,'Color', c(i,:));
                        y2 = yline(this.hAx,this.data.constraintFunctions(i).lbound,'Color', c(i,:));


                        set(get(get(y1,'Annotation'),'LegendInformation'), 'IconDisplayStyle', 'off');
                        set(get(get(y2,'Annotation'),'LegendInformation'), 'IconDisplayStyle', 'off');

                    end
                end

                legend(this.hAx, [this.leg_constr]);

                drawnow;
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