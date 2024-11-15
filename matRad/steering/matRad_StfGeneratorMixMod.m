classdef matRad_StfGeneratorMixMod < matRad_StfGeneratorBase

    properties (Constant)
        name = 'Mix Modality stf Generator';
        shortName = 'MixModStf';
        possibleRadiationModes = {'MixMod'};
    end 

    properties
        singleModalityStfs;
    end
    
    methods 
        function this = matRad_StfGeneratorMixMod(pln)
            if nargin < 1
                pln = [];
            end
            this@matRad_StfGeneratorBase(pln);

            if isempty(this.radiationMode)
                this.radiationMode = 'MixMod';
            end

            for modalityIdx=1:pln.numOfModalities
                singleModStf = matRad_StfGeneratorBase.getGeneratorFromPln(pln.originalPlans(modalityIdx));
                this.singleModalityStfs = [this.singleModalityStfs, singleModStf];
            end
        end            
    end

    methods (Access = protected)        
        function pbMargin = getPbMargin(this)
            pbMargin = min([this.singleModalityStfs(:).bixelWidth]);
        end        
    end

    methods (Static)
        function [available,msg] = isAvailable(pln,machine)
    
            %checkBasic
            try    
                %check modality
                checkModality = any(any(strcmp(matRad_StfGeneratorMixMod.possibleRadiationModes, pln.radiationMode)));
     
                preCheck = checkModality;
                
                if ~preCheck
                    return;
                end
            catch
                msg = 'Your machine file is invalid and does not contain the basic field (meta/data/radiationMode)!';
                return;
            end

            available = preCheck;
        end
    end
end
