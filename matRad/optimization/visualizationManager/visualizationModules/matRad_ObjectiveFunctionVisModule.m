classdef matRad_ObjectiveFunctionVisModule < matRad_VisualizationModule

    properties
        data;
        leg;
        leg_constr;

        dataRequest = {'fIndv', 'f','c'};
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
                        elseif isa(objective, 'DoseConstraints.matRad_DoseConstraint') || isa(objective, 'OmegaConstraints.matRad_VarianceConstraint')
                            
                            if isa(objective, 'DoseConstraints.matRad_DoseConstraintFromObjective')
                                costFuncStruct.name = objective.objective.name;
                                costFuncStruct.values = [];                                
  
                            else
                                costFuncStruct.name = objective.name;
                                costFuncStruct.values = [];
                            end
                            
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
            
            % this.data.iter = 0;
            this.data.iter = 1;
        end

        function updateData(this, dataInput)
            
            if isfield(dataInput, 'f') && ~isempty(this.data.objectiveFunctions(1).values) && ~isequal(this.data.objectiveFunctions(1).values(end), dataInput.fIndv(1))
                this.data.iter = this.data.iter +1;
            end
            this.updateDataObjectives(dataInput);
        end
        
        function updateDataObjectives(this, dataInput)

            if isfield(dataInput, 'fIndv')
                objFunctionValues = dataInput.fIndv;
            else
                objFunctionValues = [];
            end

            if isfield(dataInput, 'f')
                totF = dataInput.f;
            else
                totF = [];
            end

            if isfield(dataInput, 'c')
                c = dataInput.c;
            else
                c = [];
            end

            %if numel(objFunctionValues) == numel(this.data.objectiveFunctions)
            if ~isempty(objFunctionValues)
                for i=1:numel(this.data.objectiveFunctions)
                    this.data.objectiveFunctions(i).values(this.data.iter) = objFunctionValues(i);
                
                end
                if exist('totF', 'var')
                    this.data.totFValues(this.data.iter) = totF;
                else
                    this.data.totFValues(this.data.iter) = sum([this.data.objectiveFunctions(:).values(this.data.iter)]);
                end
            else
                % for i=1:numel(this.data.objectiveFunctions)
                %     if ~isempty(this.data.objectiveFunctions(i).values)
                %         this.data.objectiveFunctions(i).values(this.data.iter) = this.data.objectiveFunctions(i).values(max(1,this.data.iter-1));
                %     %else
                %         %this.data.objectiveFunctions(i).values(this.data.iter) = 0;
                % 
                %     end
                % end

                %if isfield(this.data, 'totFValues') && ~isempty(this.data.totFValues)
                %    this.data.totFValues(this.data.iter) = this.data.totFValues(max(1,this.data.iter-1));
                %else
                    %this.data.totFValues(this.data.iter) = 0;
                %end
            
            end

            if ~isempty(c)
                cCounter = 1;
                for i=1:numel(this.data.constraintFunctions)
                    this.data.constraintFunctions(i).values(this.data.iter,:) = c(cCounter:cCounter+size(this.data.constraintFunctions(i).ubound,1)-1,:);
                    cCounter = cCounter+1;
                end
            % else
            %     for i=1:numel(this.data.constraintFunctions)
            %         this.data.constraintFunctions(i).values(this.data.iter,:) = zeros(1,size(this.data.constraintFunctions(i).ubound,1));
            %     end
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
                    if ~isempty(y)
                        % hPlot = plot(this.hAx,x,y,'x--','LineWidth',0.5,'XDataSource','x','YDataSource','y', 'Color', c(i,:), 'DisplayName',this.leg{i});
                        hPlot = plot(this.hAx,x,y,'x--','LineWidth',0.5,'XDataSource','x','YDataSource','y', 'Color', c(i,:));
                    end
                end

                if isfield(this.data, 'totFValues')
                    y = this.data.totFValues;
                    
                    % hPlot = plot(this.hAx,x,y,'x--','LineWidth',0.5,'XDataSource','x','YDataSource','y', 'Color', c(i+1,:),'DisplayName', 'totalFunction');
                    hPlot = plot(this.hAx,x,y,'x--','LineWidth',0.5,'XDataSource','x','YDataSource','y', 'Color', c(i+1,:));
                end
                %Create plot handle and link to data for faster update
                for i=1:numel(this.data.constraintFunctions)
                    
                    y = this.data.constraintFunctions(i).values;
                    if ~isempty(y)
                        % hPlot = plot(this.hAx,x,y,'x--','LineWidth',0.5,'XDataSource','x','YDataSource','y', 'Color', c(numel(this.data.objectiveFunctions) + i,:), 'DisplayName', this.data.constraintFunctions(i).name);
                        hPlot = plot(this.hAx,x,y,'x--','LineWidth',0.5,'XDataSource','x','YDataSource','y', 'Color', c(numel(this.data.objectiveFunctions) + i,:));
                    end
                end

                legend(this.hAx, [this.leg, 'total Function', this.leg_constr]);
                % if isempty(this.hAx.Legend)
                %     legend(this.hAx);
                % else
                %     this.hAx.Legend.String = [];
                % end
                    
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