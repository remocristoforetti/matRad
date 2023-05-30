%load generic matRad environment and phantom
clear all;
matRad_rc;
matRad_cfg = MatRad_Config.instance();
load('BOXPHANTOM.mat');

%% Setup

%Define generic plan info
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = 0;
pln.propStf.couchAngles   = 0;
pln.propStf.bixelWidth    = 3;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

pln.multScen = matRad_multScen(ct,'nomScen');

pln.propDoseCalc.doseGrid.resolution.x = 8; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 8; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 8; % [mm]

%% Select radiation mode and machine
pln.radiationMode   = 'protons'; % 'carbon', 'helium'

pln.machine         = 'Generic';

%% Instantiate the specific Biological model for calculation
modelName = 'constRBE';
pln.bioParam = matRad_bioModel(pln.radiationMode,modelName);
%% Generate the stf file
stf = matRad_generateStf(ct,cst,pln);

%% Select generic plan info
pln.propDoseCalc.calcLET = 1;

%% Compute the dij for the MCN model
dij = matRad_calcParticleDose(ct,stf,pln,cst,0);
%dij = matRad_calcPhotonDose(ct,stf,pln,cst,0);
%% Select optimization quantity
pln.propOpt.quantityOpt = 'physicalDose';

%% Optimize the SOBP
resultGUI_MCN = matRad_fluenceOptimization(dij,cst,pln);

%% Other models
modelNames = {'WED', 'CAR'};

resultGUI_LETmodel = [];
for k=1:size(modelNames,2)
    pln.bioParam = matRad_bioModel(pln.radiationMode,modelNames{k});
    dij = matRad_calcParticleDose(ct,stf,pln,cst,0);
    resultGUI_LETmodel = [resultGUI_LETmodel,matRad_fluenceOptimization(dij,cst,pln)];
end

%% constRBE
modelName = 'constRBE';
pln.bioParam = matRad_bioModel(pln.radiationMode,modelName);

dij_constRBE = matRad_calcParticleDose(ct,stf,pln,cst,0);

resultGUI_constRBE = matRad_fluenceOptimization(dij_constRBE,cst,pln);
