classdef matRad_baseDataGeneration_airWidening < matRad_baseDataGeneration
    properties

    end


    methods
        function obj = matRad_baseDataGeneration_airWidening()
            
            obj@matRad_baseDataGeneration();

        end

        function saveParameters(obj)
            matRad_cfg = MatRad_Config.instance();
            
            variableName = ['AirWideningSimulation', date(), obj.MCparams.sourceParticle, '.mat'];

            saveDirectory = [obj.workingDir,filesep, 'baseDataParameters'];
            if ~exist(saveDirectory, 'dir')
                mkdir(saveDirectory);
            end
            
            obj.saveDir = saveDirectory;
            
            classProperties = properties(obj);
            for propIdx =1:size(classProperties,1)
                saveStr.(classProperties{propIdx}) = obj.(classProperties{propIdx});
            end
            save([saveDirectory, filesep, variableName],'saveStr');
        end

         function writeSimulationFiles(obj)
  
            templateFileSimulationParameters = fileread(fullfile(obj.MCparams.templateDir,'simulationParameters.txt'));
            templateBasicFile                = fileread(fullfile(obj.MCparams.templateDir,'proton_basic_air_widening.txt'));
            templateScorerIncluder = fileread(fullfile(obj.MCparams.templateDir,'scorer.txt'));

            for energyIdx=1:obj.energyParams.nEnergies
    
                obj.writeSimulationParameters(templateFileSimulationParameters,energyIdx);
                obj.writeRunFiles(energyIdx);
                obj.writeScorerIncluder(templateScorerIncluder, energyIdx);
            end

            obj.writeBasicFile(templateBasicFile);

            for scorerIdx = 1:obj.scorerParams.nScorers
                scorerName = obj.scorers{scorerIdx};
                templateScorerFile = fileread(fullfile(obj.MCparams.templateDir,['scorer_',scorerName,'.txt']));
                obj.writeScorers(templateScorerFile,scorerIdx);
            end
        end
    end
end