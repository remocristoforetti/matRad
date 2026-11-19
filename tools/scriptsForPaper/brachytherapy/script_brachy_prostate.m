
%% Script: HDR brachytherapy plan on the PROSTATE phantom
clear;
matRad_rc;
matRad_cfg = MatRad_Config.instance();

load('PROSTATE.mat','ct','cst');

% cst{5,3}    = 'TARGET';
% cst{5,6}{1} = struct(DoseObjectives.matRad_SquaredUnderdosing(100,15));
% cst{5,6}{2} = struct(DoseObjectives.matRad_SquaredOverdosing(100,17.5));
% cst{5,5}.Priority = 1;

% cst{6,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(100,12));
% cst{6,5}.Priority = 2;
% 
% cst{7,3}    = 'OAR';
% cst{7,6}    = [];
% 
% cst{1,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(10,7.5));
% cst{8,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(10,7.5));
% cst{9,6}{1} = struct(DoseObjectives.matRad_MeanDose(1));
for i=1:size(cst,1)
    cst{i,5}.Priority = 10;
end

cst{6,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(100,12));
cst{6,5}.Priority = 1;

cst{7,6} = {};
cst{7,5}.Priority = 2;

cst{1,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,8));
cst{1,5}.Priority = 2;

cst{8,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,8));
cst{8,5}.Priority = 2;

cst{9,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,5));
cst{9,5}.Priority = 3;

pln.radiationMode  = 'brachy';
pln.machine        = 'HDR';
pln.bioModel       = 'none';
pln.multScen       = 'nomScen';
pln.numOfFractions = 1;

pln.propStf.needle.seedDistance = 10;
pln.propStf.needle.seedsNo      = 7;
pln.propStf.bixelWidth          = 3;
pln.propStf.template.type       = 'checkerboard';
% pln.propStf.template.root       = matRad_getTemplateRoot(ct,cst);

pln.propDoseCalc.engine            = 'TG43';
pln.propDoseCalc.TG43approximation = '2D';
pln.propDoseCalc.doseGrid.resolution.x = 5;
pln.propDoseCalc.doseGrid.resolution.y = 5;
pln.propDoseCalc.doseGrid.resolution.z = 5;

pln.propOpt.optimizer = 'IPOPT';

stf       = matRad_generateStf(ct,cst,pln);
%%

dij       = matRad_calcDoseInfluence(ct,cst,stf,pln);

%%
matRad_cfg.defaults.propOpt.maxIter = 1500;
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

doseBrachy = resultGUI.physicalDose;
w          = resultGUI.w;
%%
matRad_cfg = MatRad_Config.instance();
outDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','brachytherapy');
if ~isfolder(outDir)
    mkdir(outDir);
end
% save(fullfile(outDir,'brachy_prostate.mat'),'ct','cst','pln','stf','w','doseBrachy');
load(fullfile(outDir, 'brachy_prostate.mat'));
isoIdx = matRad_world2cubeIndex(matRad_getIsoCenter(cst,ct),ct);
figure('Color','w');
matRad_plotSlice(ct,'axesHandle',gca,'cst',cst,'dose',doseBrachy,'plane',3,'slice',isoIdx(3), 'doseWindow', [0,14]);
