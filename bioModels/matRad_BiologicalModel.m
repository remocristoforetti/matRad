classdef (Abstract) matRad_BiologicalModel < handle
   
    properties
       bioOpt;           % boolean indicating biological optimization (=true) or physical optimization (=false)
       quantityOpt;      % quantity used for optimizaiton

       quantityVis;      % quantity used per default for visualization
       baseData;          % machine
    end
        
    properties (SetAccess = protected)
       radiationMode;    % radiation modality
       model;            % upper case short notation of the current model (e.g. LEM)
       RBE;              % constant RBE
    end

    properties (Hidden)
       identifier;       % upper case short notation of the current model in combination with the quantity used for optimization (e.g. LEM_RBExD) probably not needed
       description;      % short description of the biological model
       AvailableQuantitiesForVis = {'physicalDose', 'RBExD'};
       AvailableModels           = {'none','constRBE','MCN','WED','CAR','LEM','HEL'};   % cell array determines available models - if cell is deleted then the corersponding model can not be generated
       state;
       PropertyToChange;
       input;
    end
       
    methods
            %Constructor
            function obj = matRad_BiologicalModel()
               obj.state = 1;
            end
            
            function obj = updateState(obj)
               switch obj.state
                  case 1
                     %Do nothing
                  case 2
                     obj = obj.changePropertyValue();
                  case 3
                     obj.baseData = obj.input; %maybe just set change property value
                  case 4
                     obj.bioOpt = obj.input; %maybe just set change property value
               end
               
            end
            

            function obj = changePropertyValue(obj)
               switch obj.PropertyToChange
                  case 'radiationMode'
                     obj.radiationMode = obj.input;
                  case 'quantityOpt'
                     obj.quantityOpt = obj.input;
                  case 'baseData'
                     obj.baseData = obj.input;
                  case 'bioOpt'
                     obj.bioOpt = obj.input;
               end
            end

        %%Setters        
            function set.quantityOpt(obj,value)
               matRad_cfg =  MatRad_Config.instance();
               
               switch obj.state

                  case 1
                     if ~isempty(value) && ischar(value)
                        obj.PropertyToChange = 'quantityOpt';
                        obj.input = value;
                        obj.state = 2;
                     else
                       matRad_cfg.dispError('Unable to set %s as quantity for optimization', value);
                       obj.state = 1;
                     end
                     
                  case 2
                     if any(strcmp(value, obj.AvailableQuantitiesForOpt))
                        obj.quantityOpt = value;
                        obj.quantityVis = value;
                        if ~strcmp(value,'physicalDose')
                           obj.input = obj.baseData;
                           obj.PropertyToChange = 'baseData';
                           obj.state = 2; %try to load baseData
                        else
                           %matRad_cfg.dispWarning('physicalDose does not require baseData. \n');
                           obj.input = [];
                           obj.PropertyToChange = 'baseData';
                           obj.state = 2;
                        end
                     else
                       matRad_cfg.dispWarning('Quantity %s not available for model %s. \n', value, obj.model);
                       if isempty(obj.quantityOpt)
                           matRad_cfg.dispWarning(' Setting default.');
                           obj.state = 2;
                           obj.input = obj.AvailableQuantitiesForOpt{1};
                           obj.PropertyToChange = 'quantityOpt';
                       else
                           obj.state = 1;
                       end
                     end
               end
               obj.updateState();
             end
            
            function set.radiationMode(obj,value)
               matRad_cfg =  MatRad_Config.instance();

               switch obj.state
                  case 1
                     if ~isempty(value) && ischar(value)
                        obj.PropertyToChange = 'radiationMode';
                        obj.input = value;
                        obj.state = 2;
                     else
                       matRad_cfg.dispError('Unable to set %s as radiation modality for model %s', value, obj.model);
                       obj.state = 1;
                     end
                           
                  case 2
                     if (~isempty(obj.AvailableradiationModalities) && any(strcmp(value, obj.AvailableradiationModalities)))
                        obj.radiationMode = value;
                     else
                       matRad_cfg.dispError('Unable to set %s as radiation modality for model %s', value, obj.model);
                     end
                     obj.state = 1;
               end
               obj.updateState();
            end

             function set.baseData(obj,value)
                matRad_cfg =  MatRad_Config.instance();

                switch obj.state
                   case 1
                      if ~isempty(value) && ischar(value)
                        obj.PropertyToChange = 'baseData';
                        obj.input = value;
                        obj.state = 2;
                      else
                        matRad_cfg.dispError('Unable to set baseData source for model %s', obj.model);
                        obj.state = 1;
                      end

                   case 2
                      if ~isempty(value)
                                matRad_cfg.dispInfo('Checking validity of base data ...');
                                try                                         %Try to load the machine
                                   machine_data = load([matRad_cfg.matRadRoot filesep 'basedata' filesep strcat(obj.radiationMode, '_', value) '.mat']);
                                   obj.baseData = value;
                                   obj.state = 3;
                                catch                                       %If does not work, try generic
                                   matRad_cfg.dispWarning('Could not find any valid baseData file: %s. Trying with Generic. \n', value);
                                   try
                                         machine_data = load([matRad_cfg.matRadRoot filesep 'basedata' filesep strcat(obj.radiationMode, '_', 'Generic') '.mat']);
                                         obj.baseData = 'Generic';
                                         obj.state = 3;
                                   catch %Give up
                                           matRad_cfg.dispWarning('Could not find any Generic machine file. \n');
                                           obj.baseData = [];
                                           obj.state = 1;
                                   end
                                end
                      else %If input is empty (from internal function, need to reset baseData and bioOpt)
                         obj.baseData = [];
                         obj.PropertyToChange = 'bioOpt';
                         obj.input = 0;
                         obj.state = 4;
                         if (strcmp(obj.model,'constRBE') && ~strcmp(obj.quantityOpt, 'physicalDose'))%Specific case ~strcmp(obj.radiationMode, 'photons')
                            obj.input = 1;
                         elseif ~strcmp(obj.quantityOpt,'physicalDose')
                            matRad_cfg.dispWarning('Machine Data might be missing. QuantityOpt set to %s, but bioOpt still turned off. \n', obj.quantityOpt);
                         end
                      end
                      
                   case 3
                      if ~(strcmp(obj.quantityOpt, 'physicalDose'))
                          if ~isempty(obj.RequiredBaseData)
                             machine_data = load([matRad_cfg.matRadRoot filesep 'basedata' filesep strcat(obj.radiationMode, '_', obj.baseData) '.mat']);
                             ValidbaseData = 0;
                             for k =1:size(obj.RequiredBaseData,2)
                                if strcmp(obj.RequiredBaseData{k}, 'RBEtable')
                                    if isprop(obj, 'RBEtable')
                                       %Let it go here
                                       ValidbaseData = ValidbaseData +1;
                                    end
                                elseif ~(isfield(machine_data.machine.data, obj.RequiredBaseData{k}))
                                   matRad_cfg.dispWarning('Could not find the following machine data: %s. \n',obj.RequiredBaseData{k});
                                   ValidbaseData = 0;
                                else
                                   ValidbaseData =  ValidbaseData + 1;
                                end
                             end

                             if ValidbaseData ~= size(obj.RequiredBaseData,2)

                                matRad_cfg.dispWarning('Insufficient base data provided for this configuration. Biological optimization turned down \n');
                                obj.baseData = [];
                                obj.PropertyToChange = 'bioOpt';
                                obj.input = 0;
                                obj.state = 2;
                             else
                                if isprop(obj, 'RBEtable') && isempty(obj.RBEtable)
                                   matRad_cfg.dispInfo('done \n');
                                   obj.state = 4;
                                   obj.input = 0;
                                   matRad_cfg.dispWarning('RBEtable might be missing');
                                   
                                else
                                   obj.state = 4;
                                   obj.input = 1;
                                   matRad_cfg.dispInfo('done \n');
                                end
                             end
                          else
                             matRad_cfg.dispInfo('No base data required for model %s. \n', obj.model);
                             obj.state = 2;
                             obj.PropertyToChange = 'bioOpt';
                             obj.input = 1;
                          end
                      else
                          matRad_cfg.dispInfo('BaseData not required for optimization quantity: physicalDose. \n');
                          obj.state = 2;
                          obj.PropertyToChange = 'bioOpt';
                          obj.input = 0;
                          obj.baseData = [];
                      end
                end
                obj.updateState();
             end
      
            function set.bioOpt(obj,value)
                matRad_cfg =  MatRad_Config.instance();               
               
                switch obj.state
                   case 1
                      if islogical(value)
                         obj.state = 2;
                         obj.PropertyToChange = 'bioOpt';
                         obj.input = value;
                      else
                         matRad_cfg.dispError('Parameter bioOpt should be logical');
                         obj.state = 1;
                      end
                
                   case 2
                      if logical(value)
                         if ~isempty(obj.baseData)
                            if (isprop(obj, 'RBEtable') && isempty(obj.RBEtable))
                               matRad_cfg.dispWarning('RBEtable might be missing');
                               obj.state = 4;
                               obj.input = false;
                            else
                              obj.state = 4;
                              obj.input = true;
                            end
                         else
                            if strcmp(obj.quantityOpt, 'physicalDose')
                               matRad_cfg.dispWarning('Cannot set bioOpt with physicalDose.');
                            else
                               matRad_cfg.dispWarning('Cannot set bioOpt without baseData.');
                            end
                            obj.state = 1;
                            obj.bioOpt = false;
                         end
                      else
                         matRad_cfg.dispInfo('bioOpt set off \n');
                         obj.PropertyToChange = 'quantityOpt';
                         obj.input = 'physicalDose';
                         obj.state = 2;
                      end
                 
                   case 4
                      if logical(value)
                        obj.bioOpt = true;
                        matRad_cfg.dispInfo('All parameters are correct, biological optimization for model: %s successfully set. \n', obj.model);                         
%                         obj.state = 1;
                      else
                        obj.bioOpt = false;
                      end
                      obj.state = 1;
                end
                obj.updateState();
            end


            function obj = set.quantityVis(obj,value)
                if ~isempty(value)
                   if any(strcmp(obj.quantityOpt, [{'RBExD'}, {'effect'}]))
                      obj.quantityVis = 'RBExD';
                   else
                      obj.quantityVis = obj.quantityOpt;
                   end
                else
                    obj.quantityVis = 'physicalDose';
                end
            end
            
            
            
            function obj = set.RBE(obj,value)
                obj.RBE = value;
            end
            
%             function plotExample(obj,energy, toPlot)
%                %Load data
%                matRad_cfg = MatRad_Config.instance();
%                load([matRad_cfg.matRadRoot, filesep, 'basedata', filesep, strcat(obj.radiationMode, '_', 'Generic', '.mat')]);
%                eidx = find([machine.data(:).energy]<energy, 1, 'last');
%                
%                %Retrive infos
%                data = machine.data(eidx);
%                lege = [];
%                
%                fakeCst{1} = 0;
%                fakeCst{2} = 'Fake';
%                fakeCst{3} = 'TARGET';
%                fakeCst{4} = 0;
%                for k=1:size(obj.ModelParameters,2)
%                   fakeCst{5}.(obj.ModelParameters{k}) 
%                end
%                fakeCst{5}.alphaX = default_AlphaX;
%                fakeCst{5}.betaX = default_BetaX;
%                
%                [alphaBixel,betaBixel] = obj.calcLQParameters(data.depths,data)
%                
%                %PhysicalDose
%                figure;
%                plot(data.depths,data.Z,'.-');
%                lege = [lege; {'phisicalDose'}];
%                
%                if size(toPlot,2) > 0
%                   for k=1:size(toPlot,2)
% 
%                      switch toPlot{k}
%                         case 'effect'
%                            y = 
%  
% 
%                      end
%                   end
%                end
%                xlabel('Depth [mm]');
%                ylabel('[Gy]');
%                
%             end
    end
end