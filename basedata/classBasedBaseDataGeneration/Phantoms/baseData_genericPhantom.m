classdef baseData_genericPhantom < handle
    properties
        name;
    end


    properties (SetAccess = protected)
        scorers;
    end

    methods
        function obj = baseData_genericPhantom()
            obj.name = 'GenericPhantom';
        end

        function addScorer(obj, scorer)
            matRad_cfg = MatRad_Config.instance();

            currentScorerNames = arrayfun(@(scorerIdx) obj.scorers{scorerIdx}.name, [1:length(obj.scorers)], 'UniformOutput', false);

            if ischar(scorer)

                if ~any(strcmp(currentScorerNames, scorer))
                
                    switch scorer
                        case 'DoseToMedium'
                            obj.scorers = [obj.scorers, {baseData_scorer_DoseToMedium.instance()}];
                        otherwise
                            matRad_cfg.dispError(['Scorer: ', scorer, 'unknown']);
                    end
                
                else
                    matRad_cfg.dispWarning('Scorer already present for this phantom');                   
                end

             elseif isa(scorer, 'baseData_genericPhantom')
                 if ~any(strcmp(currentScorerNames, scorer.name))
                    obj.scorers = [obj.scorers, scorer];
                 else
                    matRad_cfg.dispWarning('Scorer already present for this phantom'); 
                 end
             else
                    matRad_cfg.dispError('Cannot recognize scorer');
             end

        end

        function removeScorer(obj,scorerName)
            matRad_cfg = MatRad_Config.instance();

            if ~isempty(obj.scorers)
                currentScorerNames = arrayfun(@(scorerIdx) obj.scorers{scorerIdx}.name, 1:length(obj.scorers), 'UniformOutput', false);
                if any(strcmp(scorerName, currentScorerNames))
                    scorerIdx = find(strcmp(scorerName, currentScorerNames));
                    obj.scorers(scorerIdx) = [];
                else
                    matRad_cfg.dispError(['No scorer named: ', scorerName, ' found']);
                end
            else
                matRad_cfg.dispError('No scorers to be removed');
            end
        end

        %% Helpers
        function writePhantomParameters(~)
        
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('This function needs to be implemented by the subclass!');

        end

         function writeScorers(~)
        
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('This function needs to be implemented by the subclass!');

         end
    end
end