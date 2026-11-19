%% Script: Hong dose calculation and SOBP optimization on the SOBP phantom
clear;
matRad_rc;
matRad_cfg = MatRad_Config.instance();
matRad_cfg.defaults.propOpt.maxIter = 10000;


load(fullfile(matRad_cfg.primaryUserFolder,'dataPaper','phantomBuilder','phantom_SOBP_2.mat'),'ct','cst');

cst{1,6}{1} = DoseObjectives.matRad_SquaredDeviation(500,60);
cst{2,6}{1} = DoseObjectives.matRad_SquaredOverdosing(10,60);

pln.radiationMode   = 'protons';
pln.machine         = 'Generic';
pln.bioModel        = 'none';
pln.multScen        = 'nomScen';
pln.numOfFractions  = 30;

pln.propStf.gantryAngles = 270;
pln.propStf.couchAngles  = 0;
pln.propStf.bixelWidth   = 3;
pln.propStf.isoCenter    = matRad_getIsoCenter(cst,ct,0);

pln.propDoseCalc.engine = 'HongPB';
pln.propDoseCalc.doseGrid.resolution.x = 1;
pln.propDoseCalc.doseGrid.resolution.y = 1;
pln.propDoseCalc.doseGrid.resolution.z = 1;

pln.propOpt.runDAO = false;

stf = matRad_generateStf(ct,cst,pln);
dij = matRad_calcDoseInfluence(ct,cst,stf,pln);
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

doseSOBP = resultGUI.physicalDose;
w        = resultGUI.w;

outDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','doseCalculation','analyticalEngines');
if ~isfolder(outDir)
    mkdir(outDir);
end
save(fullfile(outDir,'protons_SOBP_analytical_2.mat'),'ct','cst','pln','stf','w','doseSOBP');
%% 
matRad_cfg = MatRad_Config.instance();

outDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','doseCalculation','analyticalEngines');
load(fullfile(outDir, 'protons_SOBP_analytical_2.mat'));
%%
cProfIdx = 25;% + [-5:5];
figure;
plot(ct.y, squeeze(mean(doseSOBP(:,150,cProfIdx), [2,3])));
grid on;
%%
matRad_plotSlice(ct, 'cst', cst, 'dose', doseSOBP);