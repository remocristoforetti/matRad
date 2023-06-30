classdef matRad_baseDataAnalysis < matRad_baseDataGeneration
    properties
        
    end


    methods
        function obj = matRad_baseDataAnalysis(mainSimulationClass)
            obj@matRad_baseDataGeneration();

            if ~isempty(mainSimulationClass)
                propertiesMainClass = properties(mainSimulationClass);
                for propertyIdx=1:size(propertiesMainClass,1)
                    obj.(propertiesMainClass{propertyIdx}) = mainSimulationClass.(propertiesMainClass{propertyIdx});
                end
            end
        end

        function performAnalysis(obj)
            for energyIdx = 1:obj.energyParams.nEnergies
                obj.loadData(energyIdx);
                obj.analysis();
                
            end
        end
        
        function loadData(obj,energyIdx) % This might go directly into the subclass

            for scorerIdx = 1:obj.scorerParams.nScorers
                switch obj.scorerParams.scorers{scorerIdx}

                    case 'PhaseSpace'
                       obj.loadPhaseSpace(energyIdx);
                    otherwise
                        fprintf('Scorer not yet implemented');
                end
            end
        end

        function saveOutput(obj)
            for variableIdx = 1:length(obj.saveVariables)
                prop = obj.saveVariables(variableIdx);
                saveStr.(prop{1}) = obj.(prop{1});
            end
            saveName = [obj.saveNamePrefix, '_', date(),'_', obj.MCparams.sourceParticle, '.mat'];
            save(saveName, 'saveStr');
        end
    end
end