classdef matRad_baseDataAnalysis < handle %matRad_baseDataGeneration
    properties
        
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
            
            for energyIdx = 1:obj.energyParams.nEnergies
                for scorerIdx = 1:obj.scorerParams.nScorers
                    data = obj.loadData(energyIdx, scorerIdx);
                    obj.analysis(energyIdx,data);
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

        function sigma = fitSingleGaussian(obj,x,y,binEdges,start,lb,ub)
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
        end





        function [sigma1,sigma2,w_out] = fitDoubleGaussian(obj,x,y,binEdges,start, lb, ub, w)
            %check that they are all columns
            if ~iscolumn(x), x = x'; end

            if ~iscolumn(y), y=y'; end
              
            if ~iscolumn(binEdges), binEdges=binEdges'; end

            
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
            if nargin<7
                w_out = fitObject.w;
            else
                w_out = w;
            end

            visBool = 0;
            if visBool
                figure;
                plot(x,y./sum(y.*dx), '.-');
                hold on;
                plot(x, fitObject(x));
                grid on;

                grid minor;
            end
        end

    end
end