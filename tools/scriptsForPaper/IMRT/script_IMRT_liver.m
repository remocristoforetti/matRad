%% Script: photon IMRT plan on the LIVER phantom
% Same beam configuration as script_VHEE_liver.m, for comparison
clear;
matRad_rc;
matRad_cfg = MatRad_Config.instance();

load('LIVER.mat','ct','cst');

pln.radiationMode  = 'photons';
pln.machine        = 'Generic';
pln.bioModel       = 'none';
pln.multScen       = 'nomScen';
pln.numOfFractions = 30;

pln.propStf.bixelWidth   = 5;
pln.propStf.gantryAngles = [35 110 180 250 325];
pln.propStf.couchAngles  = zeros(size(pln.propStf.gantryAngles));
pln.propStf.isoCenter    = matRad_getIsoCenter(cst,ct,0);

pln.propDoseCalc.engine = 'SVDPB';
pln.propDoseCalc.doseGrid.resolution.x = 3;
pln.propDoseCalc.doseGrid.resolution.y = 3;
pln.propDoseCalc.doseGrid.resolution.z = 3;

pln.propOpt.quantityOpt = 'physicalDose';
pln.propOpt.optimizer   = 'IPOPT';
pln.propOpt.runDAO      = false;

stf       = matRad_generateStf(ct,cst,pln);
dij       = matRad_calcDoseInfluence(ct,cst,stf,pln);
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

doseIMRT = resultGUI.physicalDose;
w        = resultGUI.w;

outDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','IMRT');
if ~isfolder(outDir)
    mkdir(outDir);
end
save(fullfile(outDir,'IMRT_liver.mat'),'ct','cst','pln','stf','w','doseIMRT');

isoIdx = matRad_world2cubeIndex(pln.propStf.isoCenter,ct);
figure('Color','w');
matRad_plotSlice(ct,'axesHandle',gca,'cst',cst,'dose',doseIMRT,'plane',3,'slice',isoIdx(3));
