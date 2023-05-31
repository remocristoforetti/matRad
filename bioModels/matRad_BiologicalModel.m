classdef matRad_BiologicalModel < handle

    properties
        model;
        radiationMode;
        
        machine;
        calcBioParameters = 0;
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
                        matRad_cfg.dispInfo('Biological model successfuly instantiated.\n');
                        obj.radiationMode = value;
                        obj.machine = [];
                        obj.updateStatus();
                    else
                        matRad_cfg.dispError(['Selected radiation modality is not available for model: ', obj.model, '\n']);
                    end
                else
                    matRad_cfg.dispError('No available radiation modality has been set.\n');
                end
            end
        end

        % base Data
        

        function set.machine(obj,value)

            matRad_cfg = MatRad_Config.instance();

            currmachine = obj.machine;
            if  ~isempty(value)
                try
                    machine_name = [value.meta.machine];
                catch

                    machine_name = [value.meta.name]; %photon machine has meta.name instead of meta.machine
                end
            else
                machine_name = [];
            end
            
            if  ~isempty(currmachine)
                currmachine_name = [currmachine.meta.machine];
            else
                currmachine_name = [];
            end

            if ~isempty(obj.radiationMode)
                 currmachine_name = [obj.radiationMode, '_', currmachine_name];
                 machine_name     = [obj.radiationMode, '_', machine_name];
            end
            if ~(isempty(currmachine) && isempty(value)) && ~strcmp(currmachine_name,machine_name)
                if ~isempty(value)
                    if ~isempty(obj.radiationMode) %here if radiationMode is empty, dose nothing -> change later
                                
                           matRad_cfg.dispInfo('Checking validity of base data ...');
                           
                           fieldNames = fieldnames(value.data(1));
                           if ~isempty(obj.RequiredBaseData) && ~isempty(fieldNames)
            
                               Validmachine = 0;
                               for k =1:size(obj.RequiredBaseData,2)
                                      
                                   if ~any(strcmp(fieldNames, obj.RequiredBaseData{k}))
                                      matRad_cfg.dispWarning('Could not find the following machine data: %s.\n',obj.RequiredBaseData{k});
                                      Validmachine = 0;
                                   else
                                      Validmachine =  Validmachine + 1;
                                   end
                               end
           
                               if Validmachine ~= size(obj.RequiredBaseData,2)            
                                   matRad_cfg.dispWarning('Insufficient base data provided for this configuration. Biological optimization turned down\n');
                                   obj.machine = [];
            
                               else
                                   matRad_cfg.dispInfo(' done\n');
                                   obj.machine = value;
                                   obj.updateStatus();
                               end
            
                           elseif isempty(obj.RequiredBaseData)
                                matRad_cfg.dispInfo(' done\n');
                                matRad_cfg.dispInfo(['No machine required for model: ', obj.model, '\n']);
                                obj.machine = [];
                                obj.updateStatus();
                           end
                    end
                 else % if isempty(value)
                        %The input value is empty, just reset the machine to empty
                       obj.machine = [];
                       obj.updateStatus();
                end % end ~isempty(value)
            end
       end

       function updateStatus(obj)
            
            matRad_cfg = MatRad_Config.instance();


            currmachine = obj.machine;
            obj.machine = currmachine;


            currRadiationMode = obj.radiationMode;
            currmachine = obj.machine;

            if ~isempty(currRadiationMode) && ~isempty(currmachine)

                if ~strcmp(obj.radiationMode, 'photons')
                    obj.calcBioParameters = 1;
                    matRad_cfg.dispInfo('All parameters are correct, enabling biological model calculation.\n');
                end
            else
                obj.calcBioParameters = 0;
            end

        end


    end

end