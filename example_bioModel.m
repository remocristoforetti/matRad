%load generic matRad environment and phantom
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
pln.radiationMode   = 'carbon'; % 'carbon', 'helium'

pln.machine         = 'Generic';

%% Instantiate the specific Biological model for calculation
modelName = 'LEM';
pln.bioParam = matRad_bioModel(pln.radiationMode,modelName);
%% Generate the stf file
stf = matRad_generateStf(ct,cst,pln);

%% Select generic plan info
pln.propDoseCalc.calcLET = 1;
quantityOpt = 'RBExD';

%% Compute example of LET based bioModel
%Set the biological parameters;
% modelName = 'MCN';
% pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName,pln.machine);

%% Compute the dij for the MCN model
dij_MCN = matRad_calcParticleDose(ct,stf,pln,cst,0);

%% Optimize the SOBP
resultGUI_MCN = matRad_fluenceOptimization(dij_MCN,cst,pln);