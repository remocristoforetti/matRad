matRad_rc;
matRad_cfg = MatRad_Config.instance();
load('BOXPHANTOM.mat');
%% Setup

pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = 0;
pln.propStf.couchAngles   = 0;
pln.propStf.bixelWidth    = 3;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

pln.radiationMode   = 'carbon';

pln.machine         = 'Generic';

pln.multScen = matRad_multScen(ct,'nomScen');

pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]
stf = matRad_generateStf(ct,cst,pln);

pln.propDoseCalc.calcLET = 1;

load([pln.radiationMode, '_',pln.machine]);
quantityOpt = 'RBExD';

modelName = 'tabulatedRBEModel';
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName,pln.machine);
pln.bioParam.RBEtable = 'RBEtable_rapidMKM_Kase2008_corrected_beta_Survival_rn_0.35';
pln.bioParam.weightMode = 'eD';
%% DoseCalc

dij_MKM = matRad_calcParticleDose(ct,stf,pln,cst,0);
resultGUI = matRad_fluenceOptimization(dij_MKM,cst,pln);
%% Recalculation LEM
%tabulated LEM
modelName = 'tabulatedRBEModel';
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName,pln.machine);
pln.bioParam.weightMode = 'LET';
pln.bioParam.RBEtable = 'RBEtable_rapidLEM_Russo2011_Survival_Dt_60';
dij_LEM_tab = matRad_calcParticleDose(ct,stf,pln,cst,0);

resultGUI_LEM = matRad_calcCubes(resultGUI.w,dij_LEM_tab,1);

%% Recalculate MCN
modelName = 'MKMLET';
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName,pln.machine);
dij_MCN = matRad_calcParticleDose(ct,stf,pln,cst,0);

resultGUI_MCN = matRad_calcCubes(resultGUI.w,dij_MCN,1);
%% Plots

sliceNum = 80;
profileNum = 80;

x = [1:160]*ct.resolution.y - ct.resolution.y;
MKMprofile = resultGUI.RBExD(:,profileNum,sliceNum);
LEMprofile = resultGUI_LEM.RBExD(:,profileNum,sliceNum);
MCNprofile = resultGUI_MCN.RBExD(:,profileNum, sliceNum);

figure;

subplot(1,2,1);
imagesc(resultGUI.RBExD(:,:,sliceNum));
xline(profileNum);
subplot(1,2,2);
plot(x,MKMprofile, '.-');
hold on;
plot(x,LEMprofile,'.-');
plot(x,MCNprofile,'.-');
grid on;

grid minor;
legend('MKM', 'LEM', 'MCN');
