%% Script: RBExDose optimization with the McNamara model on the TG166 brain phantom
clear;
matRad_rc;
matRad_cfg = MatRad_Config.instance();


load(fullfile(matRad_cfg.primaryUserFolder,'patients','TG166_brain.mat'),'ct','cst');

%% Objectives (RBExDose, total over all fractions)
% the phantom only ships with target objectives, set all objectives here
cst(:,6) = {[]};
getStructIdx = @(name) find(strcmp(cst(:,2),name));

prescription = 60; % [Gy(RBE)]

cst{getStructIdx('PTV'),6}{1}         = DoseObjectives.matRad_SquaredDeviation(800,prescription);
cst{getStructIdx('GTV'),6}{1}         = DoseObjectives.matRad_SquaredDeviation(800,prescription);
cst{getStructIdx('Brain Stem'),6}{1}  = DoseObjectives.matRad_SquaredOverdosing(300,50);
cst{getStructIdx('OpticChiasm'),6}{1} = DoseObjectives.matRad_SquaredOverdosing(300,50);
cst{getStructIdx('R OpticNerve'),6}{1} = DoseObjectives.matRad_SquaredOverdosing(300,50);
cst{getStructIdx('L OpticNerve'),6}{1} = DoseObjectives.matRad_SquaredOverdosing(300,50);
cst{getStructIdx('R Eye'),6}{1}       = DoseObjectives.matRad_SquaredOverdosing(100,30);
cst{getStructIdx('L Eye'),6}{1}       = DoseObjectives.matRad_SquaredOverdosing(100,30);
cst{getStructIdx('patient'),6}{1}     = DoseObjectives.matRad_SquaredOverdosing(10,40);

%% Plan
pln.radiationMode   = 'protons';
pln.machine         = 'Generic';
pln.bioModel        = 'MCN';
pln.multScen        = 'nomScen';
pln.numOfFractions  = 30;

pln.numOfFractions       = 22;
pln.propStf.gantryAngles = [270 225 35];
pln.propStf.couchAngles   = [0 0 90];
pln.propStf.bixelWidth   = 3;
pln.propStf.isoCenter    = matRad_getIsoCenter(cst,ct,0);

pln.propDoseCalc.engine = 'HongPB';
pln.propDoseCalc.doseGrid.resolution.x = 3;
pln.propDoseCalc.doseGrid.resolution.y = 3;
pln.propDoseCalc.doseGrid.resolution.z = 3;

pln.propOpt.quantityOpt = 'RBExDose';
pln.propOpt.optimizer   = 'IPOPT';
pln.propOpt.runDAO      = false;

%% Dose calculation and optimization
stf       = matRad_generateStf(ct,cst,pln);
dij       = matRad_calcDoseInfluence(ct,cst,stf,pln);
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

physicalDose = resultGUI.physicalDose;
RBExDose     = resultGUI.RBExDose;
LET          = resultGUI.LET;
w            = resultGUI.w;

outDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','bioModels');
if ~isfolder(outDir)
    mkdir(outDir);
end
save(fullfile(outDir,'protons_TG166brain_MCN.mat'),'ct','cst','pln','stf','w', ...
    'physicalDose','RBExDose','LET');

%% Plots
isoIdx = matRad_world2cubeIndex(pln.propStf.isoCenter,ct); % [yIdx xIdx zIdx]

figure('Color','w');
subplot(1,2,1);
matRad_plotSlice(ct,'axesHandle',gca,'cst',cst,'dose',RBExDose,'plane',3,'slice',isoIdx(3), ...
    'doseWindow',[0 1.1*prescription/pln.numOfFractions],'title','RBExDose MCN [Gy(RBE)]');
subplot(1,2,2);
matRad_plotSlice(ct,'axesHandle',gca,'cst',cst,'dose',LET,'plane',3,'slice',isoIdx(3), ...
    'title','LET_d [keV/\mum]');

resultGUI = matRad_planAnalysis(resultGUI,ct,cst,stf,pln);
