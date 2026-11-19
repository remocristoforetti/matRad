%% Script: VHEE plan on the LIVER phantom
clear;
matRad_rc;
matRad_cfg = MatRad_Config.instance();

load('LIVER.mat','ct','cst');

pln.radiationMode  = 'VHEE';
pln.machine        = 'Generic';
pln.bioModel       = 'none';
pln.multScen       = 'nomScen';
pln.numOfFractions = 30;

pln.propStf.energy       = 200;
pln.propStf.bixelWidth   = 5;
pln.propStf.gantryAngles = [35 110 180 250 325];
pln.propStf.couchAngles  = zeros(size(pln.propStf.gantryAngles));
pln.propStf.isoCenter    = matRad_getIsoCenter(cst,ct,0);

pln.propDoseCalc.engine = 'HongPB';
pln.propDoseCalc.doseGrid.resolution.x = 3;
pln.propDoseCalc.doseGrid.resolution.y = 3;
pln.propDoseCalc.doseGrid.resolution.z = 3;

pln.propOpt.quantityOpt = 'physicalDose';

pln.propOpt.optimizer   = 'IPOPT';
pln.propOpt.runDAO      = false;
stf       = matRad_generateStf(ct,cst,pln);
dij       = matRad_calcDoseInfluence(ct,cst,stf,pln);
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

doseVHEE = resultGUI.physicalDose;
w        = resultGUI.w;

outDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','VHEE');
if ~isfolder(outDir)
    mkdir(outDir);
end
save(fullfile(outDir,'VHEE_liver.mat'),'ct','cst','pln','stf','w','doseVHEE');

isoIdx = matRad_world2cubeIndex(pln.propStf.isoCenter,ct);
figure('Color','w');
matRad_plotSlice(ct,'axesHandle',gca,'cst',cst,'dose',doseVHEE,'plane',3,'slice',isoIdx(3));
