classdef matRad_LexicographicEngine < handle

    properties

        iteration;
        list;
        pln;
        cst;
        dij;
        saveIntermediateSteps;
        wInitFirst;

    end

    properties (SetAccess = protected)
        currentList;
        optimizedList;
        finalWeights;
    end

    properties
        savePath;
        saveFileName = 'lex_';
    end

    methods

        function this = matRad_LexicographicEngine(iteration, list, pln, cst, dij)

            % Validate the inputs later
            this.iteration = iteration;
            this.list      = list.createListCopy();
            this.pln       = pln;
            this.cst       = cst;
            this.dij       = dij;

            %this.currentList = copy(list); % List at current step of optimization
            matRad_cfg = MatRad_Config.instance();

            if isfield(pln, 'propOpt') && isfield(pln.propOpt, 'saveDir')
                this.savePath = pln.propOpt.saveDir;
            else
                this.savePath = fullfile(matRad_cfg.primaryUserFolder, 'LexicographicProblem');
            end

        end

        function optimizeList(this)

            % Get the num of iterations needed for this list
            nIterations = this.iteration.getNumOfIterations(this.list);

            % Initialization weights for the first step, if not set 
            % matRad_initOptimization will generate the initial weights
            wInit = this.wInitFirst;

            % Perform all steps of the iteration
            for iter=1:nIterations
                                
                % Get the current single objective problem, for the current
                % iteration
                [cDij,cCst,cPln,wInit,optiProb,problemMeta] = this.iteration.getProblem(this.list, this.dij, this.cst, this.pln, wInit);

                if isempty(problemMeta) || ~isfield(problemMeta, 'skipOptimization')
                    skipOptimization = false;
                else
                    skipOptimization = problemMeta.skipOptimization;
                end

    
                % Optimize the single iteration
                [wOpt, optimizedValues, info] = this.optimizeProblem(cDij,cCst,cPln,wInit,optiProb,skipOptimization);

                this.list = this.iteration.updateList(this.list,optimizedValues, problemMeta);
                

                % Update initial weights
                wInit = wOpt;

                this.saveStep(wOpt, this.list, cPln, cCst, info, problemMeta, problemMeta.stepPriority);

            end

            this.optimizedList = this.list;
            this.finalWeights = wOpt;
        end


        function [wOpt, optimizedValues, info] = optimizeProblem(this,dij,cst,pln,wInit,optiProb,skipOptimization)

            if skipOptimization
                wOpt = wInit;

                objectiveValues = matRad_objectiveFunctions(optiProb,wInit,dij,cst);
                optimizedValues = objectiveValues;

                info = [];
               
            else % If goals not already met, need to run full optimization

                if ~isfield(pln.propOpt,'optimizer')
                    pln.propOpt.optimizer = 'IPOPT';
                end
                
                switch pln.propOpt.optimizer
                    case 'IPOPT'
                        optimizer = matRad_OptimizerIPOPT;
                    case 'fmincon'
                        optimizer = matRad_OptimizerFmincon;
                    otherwise
                        warning(['Optimizer ''' pln.propOpt.optimizer ''' not known! Fallback to IPOPT!']);
                        optimizer = matRad_OptimizerIPOPT;
                end
                        
                if ~optimizer.IsAvailable()
                    matRad_cfg.dispError(['Optimizer ''' pln.propOpt.optimizer ''' not available!']);
                end
                
                if isfield(pln.propOpt, 'acceptTollerance') && ~isempty(pln.propOpt.acceptTollerance)
                    optimizer.options.acceptable_obj_change_tol = pln.propOpt.acceptTollerance;
                end

                optimizer = optimizer.optimize(wInit,optiProb,dij,cst);

                wOpt = optimizer.wResult;
                info = optimizer.resultInfo;

                optimizedValues = matRad_objectiveFunctions(optiProb,wOpt,dij,cst);
            end

            
        end

        function saveStep(this, wOpt, list, pln, cst, info, problemMeta, iter)

            if ~exist(this.savePath, 'dir')
                mkdir(this.savePath);
            end

            fileName = fullfile(this.savePath, sprintf('%s%s_step%d.mat', this.saveFileName, this.iteration.name, iter));
            
            list = list.convertToStruct();

            save(fileName, 'wOpt', 'list', 'pln', 'cst', 'info', 'problemMeta', 'iter');
            
        end
    end

    methods (Static)
    
        
        function [results, list] = resultsFromStep(dij,varargin)
            % This function loads the saved step information and
            % reconstructs the cubes.
            % Step path can be a single  file or a full folder, in the
            % latter case, all steps in the folder are loaded and
            % reconstructed
    
            matRad_cfg = MatRad_Config.instance();


            if nargin == 2 && ischar(varargin{1}) % This if passing a dij and a path to the step to load
                if isfile(varargin{1})
                    step = matRad_LexicographicEngine.loadStepFromFile(varargin{1});
                end
            else % If passing directly the info
                % wOpt, list
                p = inputParser();

                addRequired('dij',  @(x) isstruct(x));
                addRequired('wOpt', @(x) isdouble(x));
                addOptional('list',[], @(x) isa(x, 'matRad_LexicographicList'));
                
                p = parse(varargin{:});

                dij       = p.Results.dij;
                step.wOpt = p.Results.wOpt;

            end


            results = matRad_calcCubes(step.wOpt,dij);

        end

        function step = loadStepFromFile(stepPath)

            % steppath is path to a single iteeratiokn step saved file
            step = load(stepPath);
  
        end

        function [resulGUI] = getResultsFromDir(dij, folderPath)
            
            filesInfolder = dir(folderPath);
            resulGUI = struct();

            for i=3:numel(filesInfolder)
                fileName = fullfile(filesInfolder(i).folder, filesInfolder(i).name);
                fieldName = filesInfolder(i).name(1:strfind(filesInfolder(i).name,'.')-1);
                try
                    currResult = matRad_LexicographicEngine.resultsFromStep(dij,fileName);
                catch
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispWarning(sprintf('Unable to load file: %s',strrep(fileName, '\', '\\')));
                    currResult = [];
                end
                resulGUI.(fieldName) = currResult;
            end
        end

    end
end