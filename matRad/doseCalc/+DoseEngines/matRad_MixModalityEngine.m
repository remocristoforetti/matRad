classdef matRad_MixModalityEngine < DoseEngines.matRad_DoseEngineBase

    properties (Constant)
        possibleRadiationModes = {'MixMod'};
        name = 'MixMod';
        shortName = 'MixMod';
    end

    properties
        radiationModalities;
        nModalities;
        singleModalityEngines;
        spatioTemp;
    end

    methods
        function this = matRad_MixModalityEngine(pln)
            this@DoseEngines.matRad_DoseEngineBase();

            this.assignSingleModalityEngines(pln);
        end

        function assignSingleModalityEngines(this, pln)
            this.nModalities = pln.numOfModalities;
            for modalityIdx=1:this.nModalities
                engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln.originalPlans(modalityIdx));

                this.singleModalityEngines{modalityIdx} = engine;
                this.radiationModalities{modalityIdx} = pln.originalPlans(modalityIdx).radiationMode;
                %this.spatioTemp = [this.spatioTemp, pln.];

            end

            this.spatioTemp = pln.propOpt.spatioTemp;

        end
    end

    methods (Access = protected)
        
        function dij = calcDose(this,ct,cst,stf)
            
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispInfo(' PROJECT MIXED MOD HAS BEEN ACTIVATED. AMDG. !! \n \n');

            stfModalities = {stf(:).radiationMode};
           
            for modalityIdx = 1:this.nModalities
                
                modalityName = this.radiationModalities{modalityIdx};
                
                currStf  = stf(strcmp(stfModalities,modalityName));               

                currEngine = this.singleModalityEngines{modalityIdx};
                
                dij.(modalityName) = currEngine.calcDose(ct,cst,currStf);
            end

            dij = this.finalizeDose(dij);
        end

        function dij = finalizeDose(this,dij)

            matRad_cfg = MatRad_Config.instance();

            % for modalityIdx = 1:this.nModalities
            %     modalityName = this.radiationModalities{modalityIdx};
            % 
            %     currEngine = this.singleModalityEngines{modalityIdx};
            % 
            %     dij.(modalityName) = currEngine.finalizeDose(dij.(modalityName));
            % end

            % I wanna collect all properties that might be identical among
            % the two engines into the main engine
            % Collect current single modality engine metadata

            % This might be unecessarily overcomplicated. To be simplified
            % later
            metaDataSubEngine = cellfun(@metaclass, this.singleModalityEngines, 'UniformOutput',false);
            commonProps = intersect({metaDataSubEngine{1}.PropertyList.Name}, {metaDataSubEngine{2}.PropertyList.Name});
            
            metaDataEngine = metaclass(this);
            [~,commonPropIdx] = intersect({metaDataEngine.PropertyList.Name},commonProps);

            commonProps = metaDataEngine.PropertyList(commonPropIdx);
            % Loop over properties, see if can be assigned to main
            % engine as well
            for engineProp = commonProps'
                if any(strcmp({engineProp.SetAccess}, {'public', 'protected'})) && ~engineProp.Constant

                    % get single mod properties
                    propSubEngine = cellfun(@(subEngine) subEngine.(engineProp.Name), this.singleModalityEngines, 'UniformOutput', false);
                    if isequal(propSubEngine{1}, propSubEngine{2})
                        this.(engineProp.Name) = propSubEngine{1};
                    else
                        this.(engineProp.Name) = propSubEngine;
                    end
                    
                end
            end

            % Do teh same for dij

             % for modalityIdx = 1:this.nModalities
            %     modalityName = this.radiationModalities{modalityIdx};
            allFieldsName = cellfun(@(x) fieldnames(dij.(x)), this.radiationModalities, 'UniformOutput',false);
            commonFieldsName = intersect(allFieldsName{1}, allFieldsName{2});%unique([vertcat(allFieldsName{:})]);
            
            % Exclude large fields
            [~,excludeFieldIdx] = intersect(commonFieldsName, {'physicalDose', 'mAlphaDose', 'mSqrtBetaDose', 'mLETDose'});
            commonFieldsName(excludeFieldIdx) = [];
            for propertyName = commonFieldsName'
                if isequal(dij.(this.radiationModalities{1}).(propertyName{1}), dij.(this.radiationModalities{2}).(propertyName{1}))
                    dij.(propertyName{1}) = dij.(this.radiationModalities{1}).(propertyName{1});
                else
                    dij.(propertyName{1}) = cellfun(@(x) dij.(x).(propertyName{1}),this.radiationModalities, 'UniformOutput', false);
                end
            end

            dij.numOfModalities = this.nModalities;
            dij.radiationModalities = this.radiationModalities;

            dij.spatioTemp = this.spatioTemp;
            dij.totalNumOfBixels = sum([dij.totalNumOfBixels{:}]);
        end

        % function dij = finalizeDose(this,dij)
        % 
        %     matRad_cfg = MatRad_Config.instance();
        % 
        %     stfModalities = {stf(:).radiationMode};
        % 
        %     for modalityIdx = 1:this.nModalities
        %         modalityName = this.radiationModalities{modalityIdx};
        % 
        %         currStf  = stf(strcmp(stfModalities,modalityName));               
        % 
        %         currEngine = this.singleModalityEngines{modalityIdx};
        % 
        %         tmpDij = currEngine.finalizeDose(dij.(mod   ));
        %     end
        %     %dij = finalizeDose@DoseEngines.matRad_DoseEngineBase(this,dij);
        % end
    end

    methods (Static)
        function [available,msg] = isAvailable(pln,machine)   
            % see superclass for information
           available = true;
        end
    end

end