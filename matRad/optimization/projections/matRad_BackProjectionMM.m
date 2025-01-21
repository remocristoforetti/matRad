classdef matRad_BackProjectionMM < handle
% matRad_BackProjection superclass for all backprojection algorithms 
% used within matRad optimzation processes
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    properties (SetAccess = protected)
        wCache
        wGradCache  %different cache for optimal performance (if multiple evaluations of objective but not gradient are required)
        wGradCacheProb
        d
        wGrad
        wGradProb
        dExp
        dOmegaV
        vTot
    end
    
    properties 
        dij                     %reference to matRad dij struct (to enable local changes)
        scenarios    = 1        %Scenario indices to evaluate (used for 4D & robust/stochastic optimization)
        scenarioProb = 1        %Probability associated with scenario (for stochastic optimization)
        nominalCtScenarios = 1; %nominal ct scenario (no shift, no range error) indices to evaluate (used for 4D & robust/stochastic optimization, when at least one cst structure does not have robustness)
        nModalities = 1;
        radiationModalities;
        spatioTemporalFractions;
        %totalNumOfFractions;
        useStructsForOmega = [];
    end

    
    methods
        function obj = matRad_BackProjectionMM()
            obj.wCache = [];
            obj.wGradCache = [];
            obj.wGradCacheProb = [];
            obj.d = [];
            obj.dExp = [];
            obj.dOmegaV = [];
            obj.wGrad = [];            
            obj.wGradProb = [];

        end       
        
        function obj = compute(obj,dij,w)
            if ~isequal(obj.wCache,w)
                obj.d = obj.computeResult(dij,w);
                [obj.dExp,obj.dOmegaV, obj.vTot] = obj.computeResultProb(dij,w);
                obj.wCache = w;
            end
        end
        
        function obj = computeGradient(obj,dij,doseGrad,w)
            if ~isequal(obj.wGradCache,w)
                obj.wGrad = obj.projectGradient(dij,doseGrad,w);
                obj.wGradCache = w;
            end
        end
        
        function obj = computeGradientProb(obj,dij,doseGrad,vOmegaGrad,w)
            if ~isequal(obj.wGradCacheProb,w)
                obj.wGradProb = obj.projectGradientProb(dij,doseGrad,vOmegaGrad,w);
                obj.wGradCacheProb = w;
            end
        end
        
        function d = GetResult(obj)
            d = obj.d;
        end
        
        function [dExp,dOmegaV, vTot] = GetResultProb(obj)
            dExp = obj.dExp;
            dOmegaV = obj.dOmegaV;
            vTot = obj.vTot;
        end

        function wGrad = GetGradient(obj)
            wGrad = obj.wGrad;
        end
        
        function wGrad = GetGradientProb(obj)
            wGrad = obj.wGradProb;
        end
        
        function d = computeResult(obj,dij,w)

            nBixels = cellfun(@(modality) dij.(modality).totalNumOfBixels, obj.radiationModalities);

            w = obj.splitWeigths(w, nBixels);

            for modalityIdx=1:obj.nModalities
                modalityName = obj.radiationModalities{modalityIdx};
                dTmp.(modalityName) = cell(size(dij.(modalityName).physicalDose));
                dTmp.(modalityName)(obj.scenarios) = arrayfun(@(scen) computeSingleScenario(obj,dij.(modalityName),scen,w.(modalityName)),obj.scenarios,'UniformOutput',false);
                dTmp.(modalityName)(obj.scenarios) = arrayfun(@(scen) dTmp.(modalityName){scen}*[obj.spatioTemporalFractions.(modalityName)]', obj.scenarios,'UniformOutput',false);
            end

            d = cell(size(dij.(modalityName).physicalDose));%repmat({zeros(dij.doseGrid.numOfVoxels,1)}, size(dij.(modalityName).physicalDose));
            d(obj.scenarios) = {zeros(dij.doseGrid.numOfVoxels,1)};

            for scenIdx=obj.scenarios

                for modalityIdx=1:obj.nModalities

                    modalityName = obj.radiationModalities{modalityIdx};
                    d{scenIdx} = d{scenIdx} + dTmp.(modalityName){scenIdx};
                end
            end

        end
        
        function [dExp,dOmegaV,vTot] = computeResultProb(obj,dij,w)
            % should not work with multiple CTphases at the same
            % time            
            nBixels = cellfun(@(modality) dij.(modality).totalNumOfBixels, obj.radiationModalities);

            w = obj.splitWeigths(w, nBixels);

            for modalityIdx=1:obj.nModalities
                modalityName = obj.radiationModalities{modalityIdx};

                if isfield(dij.(modalityName),'physicalDoseExp') || (isfield(dij.(modalityName), 'mAlphaDoseExp') && isfield(dij.(modalityName), 'mSqrtBetaDoseExp'))
                    if (isfield(dij.(modalityName),'physicalDoseExp') && ~isempty(dij.(modalityName).physicalDoseExp)) && ~(isfield(dij.(modalityName), 'mAlphaDoseExp') && isfield(dij.(modalityName), 'mSqrtBetaDoseExp'))
    
                        scensToInclude = find(~cellfun(@isempty, dij.(modalityName).physicalDoseExp));
                        
                        dTmpExp.(modalityName)    = cell(size(dij.(modalityName).physicalDoseExp));
                        dTmpOmegaV.(modalityName) = cell(size(dij.(modalityName).physicalDoseOmega));
                        vTotTmp.(modalityName)    = cell(size(dij.(modalityName).physicalDoseOmega));
                    
                    elseif (isfield(dij.(modalityName), 'mAlphaDoseExp') && ~isempty(dij.(modalityName).mAlphaDoseExp))
                        scensToInclude = find(~cellfun(@isempty, dij.(modalityName).mAlphaDoseExp));
        
                        dTmpExp.(modalityName)    = cell(size(dij.mAlphaDoseExp));
                        dTmpOmegaV.(modalityName) = cell(size(dij.mAlphaDoseOmega));
                        vTotTmp.(modalityName)    = cell(size(dij.mAlphaDoseOmega));
                    end
        
                       
                    [dTmpExp.(modalityName), ...
                    dTmpOmegaV.(modalityName),...
                    vTotTmp.(modalityName)] = arrayfun(@(scen) computeSingleScenarioProb(obj,dij.(modalityName),scen,w.(modalityName)),scensToInclude, 'UniformOutput',false);
                
                    dTmpExp.(modalityName)    = arrayfun(@(scen) dTmpExp.(modalityName){scen}    * [obj.spatioTemporalFractions.(modalityName)]',scensToInclude,'UniformOutput',false);
                    dTmpOmegaV.(modalityName){1}(obj.useStructsForOmega) = arrayfun(@(struct) dTmpOmegaV.(modalityName){1}{struct} * [obj.spatioTemporalFractions.(modalityName).^2]',obj.useStructsForOmega','UniformOutput',false);
                    vTotTmp.(modalityName){1}(obj.useStructsForOmega)    = arrayfun(@(struct) vTotTmp.(modalityName){1}{struct}    * [obj.spatioTemporalFractions.(modalityName).^2]',obj.useStructsForOmega','UniformOutput',false); 

                % else
                %     dTmpExp.(modalityName) = [];
                %     dTmpOmegaV.(modalityName) = [];
                %     vTotTmp.(modalityName) = [];
                end
            end

            if exist('dTmpExp', 'var')
                dExp = zeros(dij.doseGrid.numOfVoxels,1);
                vTot = repmat({0},size(vTotTmp.(modalityName){1},1),size(vTotTmp.(modalityName){1},2));
    
                for modalityIdx=1:obj.nModalities
    
                    modalityName = obj.radiationModalities{modalityIdx};
                    dExp = dExp + dTmpExp.(modalityName){1};
                    dOmegaV.(modalityName) = dTmpOmegaV.(modalityName){1};
                    vTot(obj.useStructsForOmega) = arrayfun(@(struct) vTot{struct} + vTotTmp.(modalityName){1}{struct}, obj.useStructsForOmega, 'UniformOutput',false);
                end
    
                dExp = {dExp};
            else
                dExp = [];
                dOmegaV = [];
                vTot = [];
            end

        end
     
        function wGrad = projectGradient(obj,dij,doseGrad,w)
            nBixels = cellfun(@(modality) dij.(modality).totalNumOfBixels, obj.radiationModalities);
 
            w = obj.splitWeigths(w, nBixels);

            for modalityIdx=1:obj.nModalities
                modalityName = obj.radiationModalities{modalityIdx};
                wGradTmp.(modalityName) = cell(size(dij.(modalityName).physicalDose));
                wGradTmp.(modalityName)(obj.scenarios) = arrayfun(@(scen) projectSingleScenarioGradient(obj,dij.(modalityName),doseGrad,scen,w.(modalityName)),obj.scenarios,'UniformOutput',false);
                wGradTmp.(modalityName)(obj.scenarios) = arrayfun(@(scen) wGradTmp.(modalityName){scen}.*[obj.spatioTemporalFractions.(modalityName)],obj.scenarios,'UniformOutput', false);
            end

            wGrad = cell(size(dij.(modalityName).physicalDose));
            for scenIdx=obj.scenarios

                for modalityIdx=1:obj.nModalities

                    modalityName = obj.radiationModalities{modalityIdx};
                    wGrad{scenIdx} = [wGrad{scenIdx}; wGradTmp.(modalityName){scenIdx}(:)];
                end
            end
        end
        
        function wGrad = projectGradientProb(obj,dij,dExpGrad,dOmegaVgrad,w)

            nBixels = cellfun(@(modality) dij.(modality).totalNumOfBixels, obj.radiationModalities);

            w = obj.splitWeigths(w, nBixels);

            for modalityIdx=1:obj.nModalities
                modalityName = obj.radiationModalities{modalityIdx};

                if (isfield(dij.(modalityName),'physicalDoseExp') && ~isempty(dij.(modalityName).physicalDoseExp)) && ~(isfield(dij.(modalityName), 'mAlphaDoseExp') && isfield(dij.(modalityName), 'mSqrtBetaDoseExp'))
                    wGradTmp.(modalityName) = cell(size(dij.(modalityName).physicalDoseExp));
                elseif (isfield(dij.(modalityName), 'mAlphaDoseExp') && ~isempty(dij.(modalityName).mAlphaDoseExp))
                    wGradTmp.(modalityName) = cell(size(dij.(modalityName).mAlphaDoseExp));
                end

                wGradTmp.(modalityName) = cell(size(dij.(modalityName).physicalDoseExp));
                dTmpExpGrad = {dExpGrad{1} * obj.spatioTemporalFractions.(modalityName)};
                wGradTmp.(modalityName) = projectSingleScenarioGradientProb(obj,dij.(modalityName),dTmpExpGrad,dOmegaVgrad.(modalityName),1);
            end

            
            wGrad = [];
            for modalityIdx=1:obj.nModalities
                modalityName = obj.radiationModalities{modalityIdx};
                wGrad = [wGrad; wGradTmp.(modalityName)(:)];
            end

            wGrad = {wGrad};
            % wGrad = cell(size(dij.physicalDose));
            % wGrad(obj.scenarios) = arrayfun(@(scen) projectSingleScenarioGradientProb(obj,dij,dExpGrad,dOmegaVgrad,scen,w),obj.scenarios,'UniformOutput',false);
        end

        function splitW = splitWeigths(this,w,bixelNumbers)
            bxIdx = 1;
            for modalityIdx=1:this.nModalities
                modalityName = this.radiationModalities{modalityIdx};

                modalityWeights = w(bxIdx:numel(this.spatioTemporalFractions.(modalityName))*bixelNumbers(modalityIdx) + bxIdx-1);
                splitW.(this.radiationModalities{modalityIdx}) = reshape(modalityWeights, bixelNumbers(modalityIdx),numel(this.spatioTemporalFractions.(modalityName)));
                bxIdx = bxIdx + numel(this.spatioTemporalFractions.(modalityName))*bixelNumbers(modalityIdx);
            end

        end
    end
   
    
    %These should be abstract methods, however Octave can't parse them. As soon 
    %as Octave is able to do this, they should be made abstract again 
    methods %(Abstract)
        function d = computeSingleScenario(obj,dij,scen,w)
            error('Function needs to be implemented');
        end
        
        function wGrad = projectSingleScenarioGradient(obj,dij,doseGrad,scen,w)
            error('Function needs to be implemented');
        end
        
        function [dExp,dOmegaV,vTot] = computeSingleScenarioProb(obj,dij,scen,w)
            %warning('');
        end
        
        function [wGrad] = projectSingleScenarioGradientProb(obj,dij,dExpGrad,dOmegaVgrad,scen,w)
            %warning('');
        end
    end
    
    methods (Static)
        function optiFunc = setBiologicalDosePrescriptions(optiFunc,alphaX,betaX)
            %Does nothing in a usual normal setting but return the original
            %optiFunc
        end
    end
end

