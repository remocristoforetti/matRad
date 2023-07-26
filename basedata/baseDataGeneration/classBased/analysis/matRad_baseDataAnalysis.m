classdef matRad_baseDataAnalysis < handle %matRad_baseDataGeneration
    properties
        outputAnalysis;
        saveVariables;
        saveNamePrefix;
    
    end


    methods
        function obj = matRad_baseDataAnalysis()%(mainSimulationClass)
            % obj@matRad_baseDataGeneration();
            % 
            % if ~isempty(mainSimulationClass)
            %     propertiesMainClass = properties(mainSimulationClass);
            %     for propertyIdx=1:size(propertiesMainClass,1)
            %         obj.(propertiesMainClass{propertyIdx}) = mainSimulationClass.(propertiesMainClass{propertyIdx});
            %     end
            % end
        end

        function performAnalysis(obj)
            obj.outputAnalysis = [];

            if obj.MCparams.previousRuns > 0
                obj.mergeResults();
            end

            for energyIdx = 1:obj.energyParams.nEnergies
                for scorerIdx = 1:obj.scorerParams.nScorers
                %scorerIdx = 1;
                    data = obj.loadData(energyIdx, scorerIdx);
                    obj.analysis(energyIdx,data, scorerIdx);
                end
            
            end
                
        end
        
        function loadData(obj,energyIdx)
            % to be defined by specific subclass
        end

        function analysis(obj)

            % to be defined by specific subclass
        end

        function saveOutput(obj)
            matRad_cfg = MatRad_Config.instance();

            for variableIdx = 1:length(obj.saveVariables)
                prop = obj.saveVariables(variableIdx);
                saveStr.(prop{1}) = obj.(prop{1});
            end
            
            saveFolder = [obj.workingDir, filesep, 'output', filesep];

            if ~exist('saveFolder','dir')
                mkdir(saveFolder);
            end

            saveName = [saveFolder,obj.saveNamePrefix, '_', date(),'_', obj.MCparams.sourceParticle, '.mat'];
            
            save(saveName, 'saveStr');
        end

        function sigma = fitSingleGaussian(obj,x,y,binEdges,start,lb,ub,visBool)
            %check that they are all columns
            if ~iscolumn(x), x = x'; end

            if ~iscolumn(y), y=y'; end
              
            if ~iscolumn(binEdges), binEdges=binEdges'; end
            
            dx = binEdges(2:end) - binEdges(1:end-1);

            gauss1       = @(p,x) (1./(sqrt(2*pi)*p)) * exp(-x.^2./(2*p.^2));
            objFunc     = fittype(gauss1);
            
            if ~exist('start', 'var')
                start = 0.1;
            end
            if ~exist('lb', 'var')
                lb = 0;
            end
            if ~exist('ub', 'var')
                ub = 100;
            end

            [fitObject1] = fit(x,y./sum(y.*dx), objFunc, 'Start', start, 'Lower', lb, 'Upper', ub);
            
            sigma = fitObject1.p;

            if ~exist('visBool', 'var')
                visBool = 0;
            end
            if visBool
                figure;
                subplot(1,2,1);
                plot(x,y./sum(y.*dx), '.-');
                hold on;
                plot(x, fitObject1(x));
                grid on;
                grid minor;

                subplot(1,2,2);
                semilogy(x,y./sum(y.*dx), '.-');
                hold on;
                semilogy(x, fitObject1(x));
                grid on;
                grid minor;
                sgtitle('Sigma');
            end
        end





        function [sigma1,sigma2,w_out] = fitDoubleGaussian(obj,x,y,binEdges,start, lb, ub, w)
            %check that they are all columns
            if ~iscolumn(x), x = x'; end

            if ~iscolumn(y), y=y'; end
              
            if ~iscolumn(binEdges), binEdges=binEdges'; end

            if ~exist('w', 'var'), w = []; end
            dx = binEdges(2:end) - binEdges(1:end-1);            
            gauss1       = @(p,x) (1./(sqrt(2*pi)*p)) * exp(-x.^2./(2*p.^2));
            
            switch nargin
                case 5
                    lb = [0, 0, 0];
                    ub = [1, 100, 100];
                case 6
                    ub = [1, 100, 100];
                case 7

                case 8

                otherwise
                    start = [0.01, 0.1, 1];
                    lb = [0, 0, 0];
                    ub = [1, 100, 100];
            end

            if ~(ub(2) - lb(2) < 10^(-3)) || (ub(3) - lb(3) < 10^(-3))


                if nargin<8
                    gauss2 =  @(w, p1, p2, x) ((1-w) * gauss1(p1,x) + w * gauss1(p2,x));
                    objFunc     = fittype(gauss2);
                else
                    gauss2    = @(p1, p2, x) ((1-w).* gauss1(p1,x) + w.* gauss1(p2,x));
                    objFunc     = fittype(gauss2,'coeff',{'p1','p2'}, 'independent', {'x'});
    
                end
                
                [fitObject] = fit(x,y./sum(y.*dx), objFunc, 'Start', start, 'Lower', lb, 'Upper', ub);
                
                sigma1 = fitObject.p1;
                sigma2 = fitObject.p2;
                % if isempty(w)
                %     w_out = fitObject.w;
                % else
                %     w_out = w;
                % end


                sigma1 = fitObject.p1;
                sigma2 = fitObject.p2;

                if sigma1 > sigma2
                    if isempty(w)
                        tmp = [sigma1, sigma2, fitObject.w];
                        sigma1 = tmp(2);
                        sigma2 = tmp(1);
                        w_out      = 1 - tmp(3);
                    end
                else
                    if isempty(w)
                        w_out = fitObject.w;
                    else
                        w_out = w;
                    end

                end

                visBool = 0;
                if visBool
                    figure;
                    subplot(1,2,1);
                    plot(x,y./sum(y.*dx), '.-');
                    hold on;
                    plot(x, fitObject(x));
                    grid on;
                    grid minor;

                    subplot(1,2,2);
                    semilogy(x,y./sum(y.*dx), '.-');
                    hold on;
                    semilogy(x, fitObject(x));
                    grid on;
                    grid minor;
                end
            else
                sigma1 = 0;
                sigma2 = 0;
                w_out = 0;
            end
        end

         function sigma = fitSingleRadialGaussian(obj,x,y,binEdges,start,lb,ub)
            %check that they are all columns
            if ~iscolumn(x), x = x'; end

            if ~iscolumn(y), y=y'; end
              
            if ~iscolumn(binEdges), binEdges=binEdges'; end
            
            dx = binEdges(2:end) - binEdges(1:end-1);

            gauss1      = @(p,x) 1./(2*pi*p.^2) * exp(-x.^2./(2*p.^2));
            objFunc     = fittype(gauss1);
            
            if ~exist('start', 'var')
                start = 0.1;
            end
            if ~exist('lb', 'var')
                lb = 0;
            end
            if ~exist('ub', 'var')
                ub = 100;
            end

            if ~(ub - lb < 10^(-3))
            
                [fitObject1] = fit(x,y, objFunc, 'Start', start, 'Lower', lb, 'Upper', ub);
                
                sigma = fitObject1.p;
    
                visBool = 0;
                if visBool
                    figure;
                    subplot(1,2,1);
                    plot(x,y, '.-');
                    hold on;
                    plot(x, fitObject1(x));
                    grid on;
                    grid minor;
    
                    subplot(1,2,2);
                    semilogy(x,y, '.-');
                    hold on;
                    semilogy(x, fitObject1(x));
                    grid on;
                    grid minor;
                    sgtitle('Sigma');
                end
            else
                sigma = start;
            end
         end

         function [sigma1,sigma2,w_out] = fitDoubleRadialGaussian(obj,x,y,binEdges,start, lb, ub, w, visBool)
            %check that they are all columns
            if ~iscolumn(x), x = x'; end

            if ~iscolumn(y), y=y'; end
              
            if ~iscolumn(binEdges), binEdges=binEdges'; end

            
            dx = binEdges(2:end) - binEdges(1:end-1);            
            gauss1       = @(p,x) 1./(2*pi*p.^2) * exp(-x.^2./(2*p.^2));
            
            switch nargin
                case 5
                    lb = [0, 0, 0];
                    ub = [1, 100, 100];
                case 6
                    ub = [1, 100, 100];
                case 7

                case 8

                case 9

                otherwise
                    start = [0.01, 0.1, 1];
                    lb = [0, 0, 0];
                    ub = [1, 100, 100];
            end

            if ~(ub(2) - lb(2) < 10^(-3)) || (ub(3) - lb(3) < 10^(-3))

                if isempty(w)
                    gauss2 =  @(w, p1, p2, x) ((1-w) * gauss1(p1,x) + w * gauss1(p2,x));
                    objFunc     = fittype(gauss2);
                else
                    gauss2    = @(p1, p2, x) ((1-w).* gauss1(p1,x) + w.* gauss1(p2,x));
                    objFunc     = fittype(gauss2,'coeff',{'p1','p2'}, 'independent', {'x'});
    
                end
                
                [fitObject] = fit(x,y, objFunc, 'Start', start, 'Lower', lb, 'Upper', ub);
                
                sigma1 = fitObject.p1;
                sigma2 = fitObject.p2;

                if sigma1 > sigma2
                    if isempty(w)
                        tmp = [sigma1, sigma2, fitObject.w];
                        sigma1 = tmp(2);
                        sigma2 = tmp(1);
                        w_out      = 1 - tmp(3);
                    end
                else
                    if isempty(w)
                        w_out = fitObject.w;
                    else
                        w_out = w;
                    end

                end

                %visBool = 0;
                if visBool
                    figure;
                    
                    subplot(1,2,1);
                    plot(x,y, '.-');
                    hold on;
                    plot(x, fitObject(x));
                    grid on;
                    grid minor;
                    subplot(1,2,2);
                    semilogy(x,y, '.-');
                    hold on;
                    semilogy(x, fitObject(x));
                    grid on;
                    grid minor;
                    sgtitle('Double Sigma');
                end
            else
                sigma1 = start(2);
                sigma2 = start(3);
                w_out = start(1);
            end
        end

        function mergeResults(obj)
            matRad_cfg = MatRad_Config.instance();
    
            if ~isfield(obj.MCparams, 'previousRunsDir') || isempty(obj.MCparams.previousRunsDir)
                matRad_cfg.dispError('Previous runs are specified but no directory is present');
            end

            % dirCounter =1;
            % 
            % mergedDirectoryNameBase = [obj.workingDir,filesep,'MergedDirectroy', obj.scorers{1}];
            % if ~exist(mergedDirectoryNameBase, 'dir')
            %     mkdir(mergedDirectoryNameBase);
            % 
            % else
            %     mergedDirectoryName = [mergedDirectoryNameBase, '_', num2str(dirCounter)];
            %     while exist(mergedDirectoryName,'dir')
            %             dirCounter = dirCounter +1;
            %             mergedDirectoryName = [mergedDirectoryNameBase, '_', num2str(dirCounter)];
            %     end
            %     mkdir(mergedDirectoryName);
            % 
            % end
            % 
            % for rundDirectoryIdx = 1:size(obj.MCparams.previousRunsDir,1)
            %     for runIdx=1:obj.MCparams.previousRuns(rundDirectoryIdx)
            %         cp
            %     end
            % end
        end
    end

end