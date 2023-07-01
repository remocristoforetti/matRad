classdef matRad_baseDataGeneration_airWidening < matRad_baseDataGeneration
    properties

    end


    methods
        function obj = matRad_baseDataGeneration_airWidening()
            
            obj@matRad_baseDataGeneration();

        end

        % function retriveMainClass(obj, fileName)
        %     matRad_cfg = MatRad_Config.instance;
        % 
        %     try
        %         load(fileName, 'saveStr');
        %     catch
        %         matRad_cfg.dispError('No file: ',fileName,' found');
        %     end
        % 
        %     if exist('saveStr', 'var')
        %         obj.simulateEnergies = saveStr.simulateEnergies;
        %         obj.energyParams.initFocus = saveStr.energyParams.initFocus;
        %         obj.MCparams = saveStr.MCparams;
        %     end
        % end

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
    end
end