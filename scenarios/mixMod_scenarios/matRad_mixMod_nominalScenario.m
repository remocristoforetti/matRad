classdef matRad_mixMod_nominalScenario < matRad_ScenarioModel

    properties (SetAccess = protected)
        name = 'mixMod_nomScen';
        origialScenarios;

        useScen;
    end

    methods
        function obj = matRad_mixMod_nominalScenario(ct, originalScenarios)
             if nargin == 0 
                superclassArgs = {};
            else
                superclassArgs = {ct};
            end            
            obj@matRad_ScenarioModel(superclassArgs{:});
            
            obj.origialScenarios = originalScenarios;

            %TODO: We could do this automatically in the superclass
            %Octave 5 has a bug there and throws an error
            obj.updateScenarios();
        end

          function scenarios = updateScenarios(obj)
            %Scenario Probability from pdf - here it is one since only one
            %scenario exist
            %TODO: In the context of an uncertainty model, we should
            %consider assigning probability according to the model, and
            %just leaving the weight 1
            obj.scenForProb = [0 0 0 0 0];
            obj.scenProb = 1;

            %Scenario weight 
            obj.scenWeight = 1;

            %set variables
            obj.totNumShiftScen = 1;
            obj.totNumRangeScen = 1;
            obj.totNumScen = obj.numOfCtScen; 
            
            %Individual shifts
            obj.relRangeShift = 0;
            obj.absRangeShift = 0;
            obj.isoShift = [0 0 0];

            obj.maxAbsRangeShift = max(obj.absRangeShift);
            obj.maxRelRangeShift = max(obj.absRangeShift);

            %Mask for scenario selection
            obj.scenMask = true(obj.numOfCtScen,obj.totNumShiftScen,obj.totNumRangeScen);
            
            %generic code
            [x{1}, x{2}, x{3}] = ind2sub(size(obj.scenMask),find(obj.scenMask));
            obj.linearMask    = cell2mat(x);
            totNumScen    = sum(obj.scenMask(:));
            
            %Get Scenario probability
            Sigma = diag([obj.shiftSD,obj.rangeAbsSD,obj.rangeRelSD./100].^2);
            d = size(Sigma,1);
            [cs,p] = chol(Sigma);
            obj.scenProb = (2*pi)^(-d/2) * exp(-0.5*sum((obj.scenForProb/cs).^2, 2)) / prod(diag(cs));
            obj.scenWeight = obj.scenProb./sum(obj.scenProb); 
            
            %Return variable
            scenarios = [0 0 0 0 0];

            if totNumScen ~= obj.totNumScen
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Check Implementation of Total Scenario computation - given %d but found %d!',obj.totNumScen,totNumScen);
                obj.totNumScen = totNumScen;
            end

            obj.useScen = 1;
          end

          function [dis] = combineScenarios(obj, ~, scenarioDistibutions, fractions)
              %input:
              % BP: backProjection
              % distributions : di for every modality, for every scenario


              % For nominal scenario, ony one scenario
              dis = {sum(scenarioDistibutions{1}{1}*fractions{1},2) + sum(scenarioDistibutions{2}{1}*fractions{2},2)};
          end

          function [g] = computeGradientOnScenarios(obj, BP, dijs, doseGradients, wt, fractions)
              %input:
              % BP: backProjection
              % dijs : original dijs for the two modalities
              % doseGradients : dose gradient for every scenario
              % weights for each modality
              % fractions for each modality


              % For nominal scenario, ony one scenario
              gt = cell(1,1);
              g = [];
              for modalityIdx=1:length(dijs)

                  BP.computeGradient(dijs{modalityIdx},doseGradients,wt{modalityIdx});

                  gt = BP.GetGradient();
                  % if ~isempty(gt)
                  g = [g; gt{1}*fractions{modalityIdx}];
                  % else
                  %     grads = [grads; zeros(dijs{modalityIdx}.totalNumOfBixels,1)];
                  % end
              end
              g = {g};
          end

    end
end