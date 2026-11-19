%% Script 3: Set up single-PB Monte Carlo comparisons (Hong vs TOPAS/FRED/MCsquare)
%% set matRad runtime configuration
clear;
matRad_rc;
matRad_cfg = MatRad_Config.instance();

%% Load the phantom built in phantomBuilder/matRad_script01_buildPhantom.m
load(fullfile(matRad_cfg.primaryUserFolder,'dataPaper','phantomBuilder','phantom_insert.mat'),'ct','cst');
mainDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','doseCalculation','MCengines', 'protonsInsert');
if ~isfolder(mainDir)
    mkdir(mainDir);
end
%% Treatment plan meta information (identical geometry to the analytical-engine comparison)
pln.radiationMode = 'protons';
pln.machine       = 'Generic';
pln.bioModel      = 'none';    % pure physical dose comparison
pln.multScen      = 'nomScen'; % single nominal scenario, no robustness

pln.propStf.gantryAngles = 270;
pln.propStf.couchAngles  = 0;

pln.propStf.bixelWidth   = 5; % [mm], irrelevant for a single beamlet, kept for a valid stf generator setup
pln.propStf.isoCenter    = matRad_getIsoCenter(cst,ct,0); % target centroid

pln.propDoseCalc.doseGrid.resolution.x = 1;
pln.propDoseCalc.doseGrid.resolution.y = 1;
pln.propDoseCalc.doseGrid.resolution.z = 1;

%% Generate steering information for a single pencil beam
stfGenerator = matRad_StfGeneratorParticleSingleBeamlet(pln);
stf = stfGenerator.generate(ct,cst);

w = 1;
numHistoriesDirect = 1e6;

%% Save the plan information for the run/read scripts
% The read-back script loads this instead of rebuilding pln/stf, so the
% geometry is guaranteed to match the one the input files are written with.
planFile = fullfile(mainDir,'planInfo.mat');
save(planFile,'ct','cst','pln','stf','w','numHistoriesDirect');

%% TOPAS: write input files, then try to read back an already-finished run
topasWorkingDir = fullfile(mainDir, 'TOPAS');
if ~isfolder(topasWorkingDir)
    mkdir(topasWorkingDir);
end

engineTopasWrite = DoseEngines.matRad_TopasMCEngine(pln);
engineTopasWrite.numHistoriesDirect = numHistoriesDirect;
engineTopasWrite.externalCalculation = 'write';
resultTopasWrite = engineTopasWrite.calcDoseForward(ct,cst,stf,w);
movefile(fullfile(resultTopasWrite.meta.TOPASworkingDir,'*'),topasWorkingDir);

%% MCsquare: write input files, then try to read back an already-finished run
mcSquareWorkingDir = fullfile(mainDir, 'MCsquare');
if ~isfolder(mcSquareWorkingDir)
    mkdir(mcSquareWorkingDir);
end

engineMCsquareWrite = DoseEngines.matRad_ParticleMCsquareEngine(pln);

engineMCsquareWrite.workingDir = mcSquareWorkingDir;
engineMCsquareWrite.numHistoriesDirect = numHistoriesDirect;
engineMCsquareWrite.externalCalculation = 'write';
resultMCsquare = engineMCsquareWrite.calcDoseForward(ct,cst,stf,w); % default: zero-dose write-mode placeholder

%% FRED: write input files, then try to read back an already-finished run
fredWorkingDir = fullfile(mainDir, 'FRED');
if ~isfolder(fredWorkingDir)
    mkdir(fredWorkingDir);
end

engineFredWrite = DoseEngines.matRad_ParticleFREDEngine(pln);
engineFredWrite.workingDir = fredWorkingDir;
engineFredWrite.numHistoriesDirect = numHistoriesDirect;
engineFredWrite.externalCalculation = 'write';
resultFred = engineFredWrite.calcDoseForward(ct,cst,stf,w);
