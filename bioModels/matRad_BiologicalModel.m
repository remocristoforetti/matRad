classdef matRad_BiologicalModel < handle

    properties
        radiationMode;
        baseData;
        calcBioParameters = 0;
        model;
    end

%% methods
    methods
        
        %Constructor
        function obj = matRad_BiologicalModel(radiationMode) 
            obj.radiationMode = radiationMode;
        end

        %%% Setters

        % radiation Mode
        function set.radiationMode(obj,value)
            
            matRad_cfg = MatRad_Config.instance();

            currRadiationMode = obj.radiationMode;

            if ~strcmp(currRadiationMode, value)
                if isprop(obj,'AvailableRadiationModalities') && ischar(value)
                    if any(strcmp(value,obj.AvailableRadiationModalities))
                        matRad_cfg.dispInfo('Biological model successfuly instantiated');
                        obj.radiationMode = value;
                        obj.baseData = [];
                        obj.updateStatus();
                    else
                        matRad_cfg.dispError(['Selected radiation modality is not available for model: ', obj.model]);
                    end
                else
                    matRad_cfg.dispError('No available radiation modality has been set.');
                end
            end
        end

        % base Data
        
        function set.baseData(obj,value)
            
            matRad_cfg = MatRad_Config.instance();

             currBaseData = obj.baseData;
 
             if ~strcmp(currBaseData,value) && ~(isempty(currBaseData) && isempty(value))
                 if ~isempty(value)
                     if ~isempty(obj.radiationMode) %here if radiationMode is empty, dose nothing -> change
                                
                            matRad_cfg.dispInfo('Checking validity of base data ...');
                            
                            try
                               machine_name = [obj.radiationMode, '_', value];
                               machine_data = load([matRad_cfg.matRadRoot filesep 'basedata' filesep machine_name '.mat']);
                               obj.baseData = value;
                            catch %If does not work, try generic
                               matRad_cfg.dispWarning('Could not find any valid baseData file: %s. Trying with Generic. \n', value);
                               try
                                   machine_name = [obj.radiationMode, '_', 'Generic'];
                                   machine_data = load([matRad_cfg.matRadRoot filesep 'basedata' filesep machine_name '.mat']);
                                   obj.baseData = 'Generic';
                               catch %Give up
                                   matRad_cfg.dispWarning('Could not find any Generic machine file. \n');
                                   obj.baseData = [];
                               end
                            end
            
                            %If loaded some machine files, see if has all the required
                            %data
            
                            if ~isempty(obj.RequiredBaseData) && ~isempty(obj.baseData)
            
                                ValidbaseData = 0;
                                for k =1:size(obj.RequiredBaseData,2)
                                       
                                    if ~(isfield(machine_data.machine.data, obj.RequiredBaseData{k}))
                                       matRad_cfg.dispWarning('Could not find the following machine data: %s. \n',obj.RequiredBaseData{k});
                                       ValidbaseData = 0;
                                    else
                                       ValidbaseData =  ValidbaseData + 1;
                                    end
                                end
            
                                if ValidbaseData ~= size(obj.RequiredBaseData,2)
            
                                    matRad_cfg.dispWarning('Insufficient base data provided for this configuration. Biological optimization turned down \n');
                                    obj.baseData = [];
            
                                else
                                    matRad_cfg.dispInfo(' done \n');
                                    obj.updateStatus();
                                end
            
                            end
                     end
                  else % if isempty(value)
                        %The input value is empty, just reset the baseData to empty
                        obj.baseData = [];
                        obj.updateStatus();
                  end % end ~isempty(value)
             end
        end

        function updateStatus(obj)
            
            matRad_cfg = MatRad_Config.instance();

            currRadiationMode = obj.radiationMode;
            obj.radiationMode = currRadiationMode;

            currBaseData = obj.baseData;
            obj.baseData = currBaseData;

            currRadiationMode = obj.radiationMode;
            currBaseData = obj.baseData;

            if ~isempty(currRadiationMode) && ~isempty(currBaseData)
                %For now, exclude alpha/beta calculation with photons
                if ~strcmp(obj.radiationMode, 'photons')
                    obj.calcBioParameters = 1;
                    matRad_cfg.dispInfo('All parameters are correct, enabling biological model calculation');
                end
            else
                obj.calcBioParameters = 0;
            end
        end


    end
end