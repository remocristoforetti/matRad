classdef matRad_BiologicalModel < handle
%  matRad_BiologicalModel
%  This is an abstract interface class to define Biological Models for use in
%  dose calculation and plan optimization.
%  Subclasses should at least implement the methods:
% 
%   calcTissueParameters()             to asses the completeness of the model-specific tissue information
%   calcLQParameter()                  to implement the specific bixel_alpha and bixel_beta calculation algorithm
% 
% All subclasses should also declare the  properties:
%
%   'AvailableRadiationModalities'         to specify the radiation modalities to which the model validity is limited
%   'RequiredBaseData'                     to check the availability of information stored in the provided machine file
%
% constructor (Abstract)
%   matRad_BiologicalModel(radiationMode)
%
% input
%   radiationMode:                 radiation modality selected for the plan
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2023 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        model;                      % name of the implemented model
        radiationMode;              % selected radiation modality
        
        machine;                    % machine data

    end

    properties (SetAccess = protected)
        calcBioParameters = 0;          % boolean to trigger the calculation of model specific biological quantitites like dose weighted alpha/beta distributions
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
                        matRad_cfg.dispError('Selected radiation modality is not available for the selected model');
                    end
                else
                    matRad_cfg.dispError('No available radiation modality has been set.');
                end
            end
        end

        % base Data
        function set.machine(obj,value)
       
        % This function sets the machine property. In order for the
        % machine to be correctly set, it has to be consistent with the
        % required base data. If the machine is missing some field
        % entry, the instantiation of this class property fails and no
        % biological calculation can be triggered for the model.
        % If this function is called but no actual change is 
       

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
                try
                    currmachine_name = [value.meta.machine];
                catch

                    currmachine_name = [value.meta.name]; %photon machine has meta.name instead of meta.machine
                end
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
                                      matRad_cfg.dispWarning('Could not find the following machine data: %s',obj.RequiredBaseData{k});
                                      Validmachine = 0;
                                   else
                                      Validmachine =  Validmachine + 1;
                                   end
                               end
           
                               if Validmachine ~= size(obj.RequiredBaseData,2)            
                                   matRad_cfg.dispError(['Insufficient base data provided for model: ', obj.model, '. Cannot perform dose calculation']);
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

            % Also set.machine calls updateStatus(), but this step then does
            % not affect the machine setting because in set.machine nothing
            % happens if the machine is the same as the stored one
            currmachine = obj.machine;
            obj.machine = currmachine;

            %Same with this step, if the radiation mode is unchanged,
            %nothing happens
            currRadiationMode = obj.radiationMode;
            currmachine = obj.machine;

            if ~isempty(currRadiationMode) && ~isempty(currmachine)

                if ~strcmp(obj.radiationMode, 'photons')
                    % If the radiation modality and the machine have been assigned to the class, then 
                    % the calculation of biological quantities is allowed. This triggers the
                    % call to calcTissueParameters() and calcLQParameter() in matRad_calcParticleDose().
                    obj.calcBioParameters = 1;
                    matRad_cfg.dispInfo('All parameters are correct, enabling biological model calculation.\n');
                end
            else
                obj.calcBioParameters = 0;
            end

       end

       function tissueParameters = calcTissueParameters(obj)
       
       % This function extracts the relevant tissue parameters (alpha/beta values) for
       % every voxel. The function needs to be implemented by the user
       % for the specific model they want to include. The
       % tissueParameter structure produced as an output is given as an
       % input to the calcLQParameter() function and can thus store the
       % required information for the actual model calculation.
       %    input
       %        cst            
       %        numVoxels      the number of Voxels in the doseGrid
       %        stf            
       %        ctScen
       %    output
       %        tissueParameters        structure containing the necessary tissue parameters. This structure is then passed as an
       %                                input to the function calcLQParameter().
       

           matRad_cfg = MatRad_Config.instance();
           tissueParameters = [];
           matRad_cfg.dispError('Function: calcTissueParameters not implemented by the bioModel subclass!');
       
       end

        function [bixelAlpha, bixelBeta] = calcLQParameter(obj)

        % This function applies the actual biological model to compute the
        % predicted alpha and beta values for every voxel in the current
        % bixel. 
        %    input
        %        vRadDepths         equivalent depths of the voxels hit by the current bixel   
        %        baseDataEntry      the machine data entry for the energy
        %                           of the current bixel, it contaions the
        %                           depth-dependent quantities to be
        %                           sampled
        %        TissueParam        structure computed by
        %                           calcTissueParameters(). Contains tissue
        %                           related information for all the voxels
        %                           on the dose Grid
        %        ix                 indexes of the voxels hit by the
        %                           current bixel, referred to the dose
        %                           grid.
        %    output
        %        bixelAlpha         alpha values for each voxel hit by the
        %                           current bixel
        %        bixelBeta          beat values for each voxel hit by the
        %                           current bixel
        %
        

            matRad_cfg = MatRad_Config.instance();
            bixelAlpha = [];
            bixelBeta  = [];
            matRad_cfg.dispError('Function: calcLQParameter not implemented by the bioModel subclass!');
       
       end


    end

end