classdef matRad_mixMod_nominalVsnominal < matRad_ScenarioModel
%This scenario modlaity combines the nominal photon scenario with all the proton scenarios and the noinal proton scenario with all the photon scenarios    
    properties (SetAccess = protected)
        name = 'mixMod_nomVsnom';
        originalScenarios;
    
        nSamples;
        useScen;
        combinationMask;
    end

    methods
        function obj = matRad_mixMod_nominalVsnominal(ct, originalScenarios)
            if nargin == 0 
                superclassArgs = {};
            else
                superclassArgs = {ct};
            end            
        obj@matRad_ScenarioModel(superclassArgs{:});
        
        obj.originalScenarios = originalScenarios;
        
        %TODO: We could do this automatically in the superclass
        %Octave 5 has a bug there and throws an error


        obj.updateScenarios();
        end

        function set.originalScenarios(obj,originalScenarios)

            for modalityIdx = 1:size(originalScenarios,2)
                nSmodality(modalityIdx) = originalScenarios{modalityIdx}.nSamples;
            end
                        
            obj.nSamples = sum(nSmodality)-1; %-1 here to avoit double counting of nominal-nominal combination    
            

            obj.originalScenarios = originalScenarios;
        end

        function scenarios = updateScenarios(obj)
            matRad_cfg = MatRad_Config.instance();
            scenarios = zeros(obj.nSamples,2);
            
            obj.scenProb    = (1./obj.nSamples)*ones(1,obj.nSamples);
            obj.scenForProb = scenarios;

            %Scenario weight
            obj.scenWeight  = ones(obj.nSamples,1)./obj.nSamples;
            
            %set variables
            if ~all(cellfun(@(scenario) getfield(scenario,'includeNominalScenario'), obj.originalScenarios))
                matRad_cfg.dispWarning('NOt all the scenarios include the nominal scenario, results will be bayased');
            end    
                
            obj.totNumShiftScen = sum(cellfun(@(scenario) getfield(scenario,'totNumShiftScen'), obj.originalScenarios))-1; %-1 to eliminate the nominal scenario accounted twice for the two modalities?
            obj.totNumRangeScen = sum(cellfun(@(scenario) getfield(scenario,'totNumRangeScen'), obj.originalScenarios))-1;            
            obj.totNumScen = sum(cellfun(@(scenario) getfield(scenario,'totNumScen'), obj.originalScenarios))-1;
            
            %Individual shifts
            obj.relRangeShift = 0; % Dummy info here
            obj.absRangeShift = 0;
            obj.isoShift = [0 0 0];
            
            obj.maxAbsRangeShift = max(obj.absRangeShift);
            obj.maxRelRangeShift = max(obj.absRangeShift);
            


            %Mask for scenario selection
            obj.scenMask = true(obj.numOfCtScen,obj.totNumShiftScen,obj.totNumRangeScen); %this is also dummy


            obj.useScen = [1:obj.nSamples];

            combMask      = zeros(obj.nSamples,size(obj.originalScenarios,2));
            
            combMask(:,1) = [1; ones(obj.originalScenarios{2}.nSamples-1,1); [2:obj.originalScenarios{1}.nSamples]'];
            combMask(:,2) = [1; [2:obj.originalScenarios{2}.nSamples]'; ones(obj.originalScenarios{1}.nSamples-1,1)];

            obj.combinationMask = combMask;
            % for modalityIdx=1:size(obj.originalScenarios,2)
            %     combMask(2:end,modalityIdx) = [[2:obj.originalScenarios{modalityIdx}.nSamples], [2:obj.originalScenarios{3-modalityIdx}.nSamples]]';
            % end
        end


        function [dis] = combineScenarios(obj, BP, scenarioDistibutions, fractions)
             %input:
              % BP: backProjection
              % distributions : di for every modality, for every scenario

              % Combine here the plans
              dis = cell(obj.nSamples,1);

              for modalityIdx=1:size(scenarioDistibutions,2)
                  modalityScenarios{modalityIdx} = find(~cellfun(@(scenPhysicalDose) isempty(scenPhysicalDose), scenarioDistibutions{modalityIdx}));
              end

              for scenIdx = 1:obj.nSamples

                  dis{scenIdx} = scenarioDistibutions{1}{modalityScenarios{1}(obj.combinationMask(scenIdx,1))}*fractions{1} + scenarioDistibutions{2}{modalityScenarios{2}(obj.combinationMask(scenIdx,2))}*fractions{2};
              end
              % for modalityIdx = 1:size(scenarioDistibutions,2)
              % 
              %     nominalModalityDistribution = scenarioDistibutions{modalityIdx}{1}*fractions{modalityIdx};
              % 
              %     distIdx = find(~cellfun(@(scenPhysicalDose) isempty(scenPhysicalDose), scenarioDistibutions{3 - modalityIdx}));% if modalityIdx=1, choose scenDist(2) else use scenDist(1)
              % 
              %     distIdx(distIdx==1) = []; % exclude the nominal-nominal scenario
              % 
              %     for distributionIdx=1:length(distIdx)
              %         currDistribution = nominalModalityDistribution + [scenarioDistibutions{3 - modalityIdx}{distIdx(distributionIdx)}].*(fractions{3 - modalityIdx});
              %         dis = [dis; {currDistribution}];
              %     end
              % end
              % 
              % nominalDistribution = scenarioDistibutions{1}{1}*fractions{1} + scenarioDistibutions{2}{1}*fractions{2};
              % dis = [{nominalDistribution}; dis];
        end

        function [g] = computeGradientOnScenarios(obj, BP, dijs, doseGradients, wt, fractions)

            g = cell(obj.nSamples,1);
            %Need to go through this way because wGradientCache is
            %protected and need to call for each modality in sequence

            nonEmptyScenFirstModality = find(~cellfun('isempty', dijs{1}.physicalDose));
            nonEmptyScenSecondModality = find(~cellfun('isempty', dijs{2}.physicalDose));
            
            for scenarioIdx =1:obj.nSamples
                BP.scenarios = nonEmptyScenFirstModality(obj.combinationMask(scenarioIdx,1)); %This corresponds to index of the dij.physicalDose
                firstModalityDoseGradient = cell(size(dijs{1}.physicalDose));
                firstModalityDoseGradient(BP.scenarios) = doseGradients(scenarioIdx);

                BP.computeGradient(dijs{1}, firstModalityDoseGradient, wt{1});
                gTmp = BP.GetGradient();
                g_firstMod = gTmp(BP.scenarios);

                BP.scenarios = nonEmptyScenSecondModality(obj.combinationMask(scenarioIdx,2)); %This corresponds to index of the dij.physicalDose
                secondModalityDoseGradient = cell(size(dijs{2}.physicalDose));
                secondModalityDoseGradient(BP.scenarios) = doseGradients(scenarioIdx);
                BP.computeGradient(dijs{2}, secondModalityDoseGradient, wt{2});
                gTmp = BP.GetGradient();

                g_secondMod = gTmp(BP.scenarios);

                g{scenarioIdx} = [g_firstMod{1}.*fractions{1}; g_secondMod{1}.*fractions{2}];
            end
%{
            %Get scenario indexes that are coupled to nominal scenario of
            %second modality

            firstModalityIdx = find(obj.combinationMask(:,2)==1);

            BP.scenarios = find(~cellfun('isempty', dijs{1}.physicalDose)); %[1:obj.originalScenarios{1}.nSamples];
            firstModalityDoseGradient = cell(size(dijs{1}.physicalDose));
            firstModalityDoseGradient(BP.scenarios) = doseGradients(firstModalityIdx);

            BP.computeGradient(dijs{1},firstModalityDoseGradient,wt{1});
            g_firstModality_complete = BP.GetGradient();

            % Second modality regular combinations
            secondModalityIdx = find(obj.combinationMask(:,1)==1);

            BP.scenarios = find(~cellfun('isempty', dijs{2}.physicalDose));
            secondModalityDoseGradient = cell(size(dijs{2}.physicalDose));
            secondModalityDoseGradient(BP.scenarios) = doseGradients(secondModalityIdx);

            BP.computeGradient(dijs{2},secondModalityDoseGradient,wt{2});
            g_secondModality_comlete = BP.GetGradient();


            BP.scenarios = 1;
            g_firstModality_nonReg = [];
            nonRegularFirstModalityIdx = secondModalityIdx(find(secondModalityIdx~=1)); %Get indexes of doseGrad that have the nominal first modality scenario (exept for nominal-nominal)
            for scenarioIdx=nonRegularFirstModalityIdx
                BP.computeGradient(dijs{1},doseGradients(scenarioIdx), wt{1});
                tmpG = BP.GetGradient();
                g_firstModality_nonReg = [g_firstModality_nonReg, tmpG(1)];
                BP.wGradCache = [];
            end

            g_secondModality_nonReg = [];
            nonRegularSecondModalityIdx = firstModalityIdx(find(firstModalityIdx~=1)); %Get indexes of doseGrad that have the nominal first modality scenario (exept for nominal-nominal)
            for scenarioIdx=nonRegularSecondModalityIdx
                BP.computeGradient(dijs{2},doseGradients(scenarioIdx), wt{2});
                tmpG = BP.GetGradient();
                g_secondModality_nonReg = [g_secondModality_nonReg, tmpG(1)];
            end

            %for scenarioIdx=1:obj.nSamples
            g(firstModalityIdx) = [g_firstModality.*fractions{1}; g_secondModality(1).*fractions{2}];

            g(secondModalityIdx~=1) = [g_firstModality(1).*fractions{1}; g_secondModality(secondModalityIdx~=1).*fractions{2}];
%}
            %end
        end
    end
end