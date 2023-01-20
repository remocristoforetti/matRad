classdef matRad_BioScenarios < matRad_ScenarioModel
%  matRad_RandomScenarios
%  Implements randomly sampled scenarios
%
% constructor
%   matRad_RandomScenarios()
%   matRad_RandomScenarios(ct)
%
% input
%   ct:                 ct cube
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2022 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties (Hidden)
        bioModelParameters;             %structure;
        AvailableBioModels = {};
    end

    %Deprecated Properties that were used
    properties (Dependent)
        numOfShiftScen;
        numOfRangeShiftScen;
    end

    properties (SetAccess=protected)
        BioModels;
    end
    
    
    methods
        function this = matRad_BioScenarios(ct)    
            if nargin == 0 
                superclassArgs = {};
            else
                superclassArgs = {ct};
            end
            
            this@matRad_ScenarioModel(superclassArgs{:});
        end

        %% Setters & Update
        function set.bioModelParameters(this,value)
            matRad_cfg = MatRad_Config.instance();
            if isstruct(value)
                if ~(sum(isfield(value,{'radiationMode', 'quantityOpt', 'baseData'})))
                   matRad_cfg.dispWarning('bioModelParameters not complete'); 
                else
                    switch value.radiationMode
                        case 'photons'
                            AvailableModels = {'constRBE'};
                            radMode = value.radiationMode;
                        case 'protons'
                            AvailableModels = {'WED', 'MCN', 'LSM', 'MKMLET', 'MKMLET_corrected', 'constRBE', 'CAR'};
                            radMode = value.radiationMode;
                        case {'helium', 'carbon'}
                            AvailableModels = {'LEM'};
                            radMode = value.radiationMode;                        
                       otherwise
                            matRad_cfg.dispWarning('Radiation mode %s not available', value.radiationMode); 
                            AvailableModels = [];
                            radMode = [];
                    end
                    this.AvailableBioModels = AvailableModels;
                    
                    if any(strcmp(value.quantityOpt,{'RBExD','effect'}))
                        opt = value.quantityOpt;
                    else
                        matRad_cfg.dispWarning('QuantityOpt not valid, setting default');
                        opt = 'RBExD';
                    end

                    if ischar(value.baseData)
                       try
                          machine_data = load([matRad_cfg.matRadRoot filesep 'basedata' filesep strcat(value.radiationMode, '_', value.baseData) '.mat']);
                          Data = value.baseData;
                       catch
                        matRad_cfg.dispWarning('Unable to load machine data %s', value.baseData);
                        Data = [];
                       end   
                    end
                    
                    str = struct('radiationMode', radMode, ...
                                 'quantityOpt', opt, ...
                                 'baseData', Data);
                    this.bioModelParameters = str;
                    
                    this.updateScenarios();
                end
            else
              matRad_cfg.dispWarning('bioModelParameters not correct');
            end
        end
       
        function set.numOfShiftScen(this,numOfShiftScen)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The property numOfShiftScen of the scenario class will soon be deprecated! Use nSamples instead');
            
            %That's for downwards compatibility
            if ~isscalar(numOfShiftScen)
                numOfShiftScen = unique(numOfShiftScen);
            end

            this.nSamples = numOfShiftScen;
        end

        function  value = get.numOfShiftScen(this)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The property numOfShiftScen of the scenario class will soon be deprecated! Use nSamples instead');            
            value = this.nSamples;
        end

        function set.numOfRangeShiftScen(this,numOfRangeShiftScen)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The property numOfRangeShiftScen of the scenario class will soon be deprecated! Use nSamples instead');           
            this.nSamples = numOfRangeShiftScen;
        end

        function  value = get.numOfRangeShiftScen(this)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The property numOfRangeShiftScen of the scenario class will soon be deprecated! Use nSamples instead');
            value = this.nSamples;
        end
        
        function BioModelInstance = setBioModel(this,BioModel)
         
             radiationMode = this.bioModelParameters.radiationMode;
             quantityOpt   = this.bioModelParameters.quantityOpt;
             machine       = this.bioModelParameters.baseData;

             matRad_cfg = MatRad_Config.instance();
             if any(strcmp(BioModel, this.AvailableBioModels))
                 BioModelInstance = {matRad_bioModel(radiationMode,quantityOpt,BioModel,machine)};
             else
                 matRad_cfg.dispError('BioModel not valid \n');
             end
             
              if ~BioModelInstance{1}.bioOpt
                 BioModelInstance{1}.bioOpt = true;
                 if ~BioModelInstance{1}.bioOpt
                     matRad_cfg.dispError('Cannnot set bioOpt for model %s \n', BioModelInstance{1}.model);
                 end
              end
        end

    end
end