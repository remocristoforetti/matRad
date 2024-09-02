function model = matRad_bioModel(sRadiationMode, sModel, sMachine)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  matRad_bioModel
%  This is a helper function to instantiate a matRad_BiologicalModel. This
%  function currently exists for downwards compatability, as the new
%  Biological Models will follow a polymorphic software architecture
%
% call
%   matRad_bioModel(sRadiationMode, sModel)
%   matRad_bioModel(sRadiationMode, sModel, sMachine)
%
%   e.g. pln.bioModel = matRad_bioModel('protons','MCN', 'Generic')
%
% input
%   sRadiationMode:     radiation modality 'photons' 'protons' 'helium' 'carbon' 'brachy'
%   
%   sModel:             string to denote which biological model is used
%                       'none': for photons, protons, carbon                'constRBE': constant RBE for photons and protons
%                       'MCN': McNamara-variable RBE model for protons      'WED': Wedenberg-variable RBE model for protons
%                       'LEM': Local Effect Model for carbon ions
%   sMachine (optional): checks the consistency of the provided machine for
%                        the default parameters of the model. Can be used
%                        to check the completeness of the base data
%                        provided
%
% output
%   model:              instance of a biological model
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

matRad_cfg = MatRad_Config.instance();

% Look for the correct inputs
p = inputParser;
isNotOldQuantityOpt = @(x) ~any(strcmp(x, {'physicalDose', 'RBExD', 'effect', 'BED'}));
addRequired(p, 'sRadiationMode', @ischar);
addRequired(p, 'sModel',isNotOldQuantityOpt);
addOptional(p, 'sMachine', [], @ischar);

p.KeepUnmatched = true;

try
    if nargin>2
        parse(p, sRadiationMode, sModel, sMachine);
    else
        parse(p, sRadiationMode, sModel);
    end

catch ME
    if matRad_cfg.isOctave
        checkString = upper('sModel');
    else

        checkString = 'sModel';
    end
    if ~isempty(strfind(ME.message, checkString)) && ~isNotOldQuantityOpt(sModel)
        matRad_cfg.dispDeprecationWarning('quantityOpt input for the biological model will be deprecated!');
        sModel = sMachine;
        sMachine = [];
    else
        matRad_cfg.dispError(ME.message);
    end
end

if ~exist('sMachine', 'var')
    sMachine = [];
end

%Check for the available models
mainFolder        = fullfile(matRad_cfg.matRadSrcRoot,'bioModels');
userDefinedFolder = fullfile(matRad_cfg.primaryUserFolder, 'bioModels');

if ~exist(userDefinedFolder,"dir")
    folders = {mainFolder};
else
    folders = {mainFolder,userDefinedFolder};
end

availableBioModelsClassList = matRad_findSubclasses('matRad_BiologicalModel', 'folders', folders , 'includeSubfolders',true);
modelInfos = matRad_identifyClassesByConstantProperties(availableBioModelsClassList,'model');
modelNames = {modelInfos.model};

if numel(unique({modelInfos.model})) ~= numel(modelInfos)
    matRad_cfg.dispError('Multiple biological models with the same name available.');
end
            
selectedModelIdx = find(strcmp(sModel, modelNames));
            
% Create first instance of the selected model
if ~isempty(selectedModelIdx)
    tmpBioParam = modelInfos(selectedModelIdx).handle();
else
    matRad_cfg.dispError('Unrecognized biological model: %s', sModel);
end

% For the time being I do not assigne the model specific parameters, they
% can be assigned by the user later

correctRadiationModality = any(strcmp(tmpBioParam.availableRadiationModalities, sRadiationMode));

% Assume base data are correct, if not the check function will fail
correctBaseData = true;

if ~isempty(sMachine)
    try
        machineStruct = struct('radiationMode', sRadiationMode, 'machine', sMachine);
        machine = matRad_loadMachine(machineStruct);
    catch
        matRad_cfg.dispError('Could not find the following machine file: %s',[sRadiationMode, '_', sMachine]);
    end

    correctBaseData = tmpBioParam.checkBioCalcConsistency(machine);
end

if correctRadiationModality && correctBaseData
    model = tmpBioParam;
elseif ~correctRadiationModality
    matRad_cfg.dispError('Incorrect radiation modality for the required biological model');
    model = [];
elseif ~correctBaseData
    matRad_cfg.dispWarning('Insufficient base data for required biological model.');
    model = tmpBioParam;
end


end % end class definition