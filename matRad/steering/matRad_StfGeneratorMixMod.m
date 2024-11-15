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
                this.singleModalityStfs = [this.singleModalityStfs, {singleModStf}];
            end
        end


        function stf = generate(this, ct, cst)
            % Generate steering information for the given ct and cst
            % This is a base class function performing the following tasks
            % 1. Checking the input
            % 2. Initializing the patient geometry (protected properties)
            % 3. Generating the source information (and thus the "stf")
            
            % Instance of MatRad_Config class
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispInfo('matRad: Generating stf struct with generator ''%s''... ',this.name);
            
            this.ct = ct;
            this.cst = cst;

            for modalityIdx=1:numel(this.singleModalityStfs)

                this.singleModalityStfs{modalityIdx}.ct = ct;
                this.singleModalityStfs{modalityIdx}.cst = cst;

            end
            cellfun(@(modality) modality.initialize(), this.singleModalityStfs);
            cellfun(@(modality) modality.createPatientGeometry(), this.singleModalityStfs);
            tmp_stf = cellfun(@(modality) modality.generateSourceGeometry(), this.singleModalityStfs, 'UniformOutput',false);
            
            fields = cellfun(@fieldnames,tmp_stf, 'UniformOutput',false);
            totalFields = unique(vertcat(fields{:}));

            for k=1:length(tmp_stf)
                isPlanField = find(~isfield(tmp_stf{k},totalFields));
                if any(isPlanField)
                    for m=isPlanField
                        tmp_stf{k} = setfield(tmp_stf{k},{1},totalFields{m},[]);
                    end
                end
            end
   
            stf = [tmp_stf{:}];
        end
    end

    methods (Access = protected)        
        function pbMargin = getPbMargin(this)
            photonIdx = find(cellfun(@(modality) isa(modality, 'matRad_StfGeneratorPhotonIMRT'), this.singleModalityStfs), 1,'first');
            if ~isempty(photonIdx)
                pbMargin = this.singleModalityStfs{photonIdx}.bixelWidth;
            end
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
