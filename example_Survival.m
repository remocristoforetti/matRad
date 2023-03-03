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

pln.machine         = 'HITgantry';
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
pln.bioParam.RBEtable = 'ExampleMKMInaniwa2010';
 

%% Override stf
stf = OverrideStf(stf,1);
%% DoseCalc

dij_MKM = matRad_calcParticleDose(ct,stf,pln,cst,0);
%% Opt
resultGUI = matRad_fluenceOptimization(dij_MKM,cst,pln);

%% Plot
Slice =80;
Profile = 80;

RBExD_Profile_MKM = resultGUI.RBExD(:,Profile,Slice);
pD_Profile_MKM = resultGUI.physicalDose(:,Profile,Slice);

x=[1:160]*ct.resolution.x - ct.resolution.x/2;
% figure;
% subplot(1,2,1);
% imagesc(resultGUI.RBExD(:,:,Slice));
% xline(Profile);
% subplot(1,2,2);
% 
% plot(x,RBExD_Profile_MKM,'.-');
% hold on;
% plot(x,pD_Profile_MKM,'.-');
% 
% grid on;

%% Recompute LEM
pln.machine = 'Generic';
modelName = 'LEM';

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName,pln.machine);
dij_LEM = matRad_calcParticleDose(ct,stf,pln,cst,0);
resultGUI_LEM = matRad_calcCubes(resultGUI.w,dij_LEM,1);
RBExD_Profile_LEM = resultGUI_LEM.RBExD(:,Profile,Slice);
%% Plot
figure;
plot(x,RBExD_Profile_MKM,'.-');
hold on;
plot(x,RBExD_Profile_LEM,'.-');
grid on;
legend('MKM', 'LEM');
%% baseData comparison
pln.machine = 'fakeGeneric';
modelName = 'tabulatedRBEModel';
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName,pln.machine);
pln.bioParam.RBEtable = 'ExampleMKMInaniwa2010';

load('carbon_HITgantry.mat');
eIdx = 31;
depthsBase = machine.data(eIdx).depths;
Z = machine.data(eIdx).Z;

matRad_calcDoseInit;
dij.TissueParameters = pln.bioParam.calcTissueParameters(cst,dij.doseGrid.numOfVoxels, stf,1);
ix = ones(size(depthsBase,1),1);%[1:size(dij.TissueParameters.tissueClass,1)];

% modelName = 'tabulatedRBEModel';
% pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName,pln.machine);
% pln.bioParam.RBEtable = 'ExampleMKMInaniwa2010';

MKMprofile = pln.bioParam.calcLQParameter(depthsBase,machine.data(eIdx),dij.TissueParameters,ix);
LEMprofile = machine.data(eIdx).alpha(:,1);
figure;
%plot(depthsBase, Z./max(Z(1)), '.-');
%hold on;
plot(depthsBase, MKMprofile, '.-');
hold on;
plot(depthsBase, LEMprofile, '.-');
grid on;
legend('MKM', 'LEM');

%% Survival class example
IntegralS = IntegralDose_Survival();

IntegralS.survivalSourcePath = '$HOME/Survival';

IntegralS.wDir = 'D:\matRad_gitHubRemo_MKM_TOPAS\MKM';
IntegralS.survivalParameterFileName = 'MKM_1';
IntegralS.calcProperties.projectName     = 'LUT_1';
IntegralS.calcProperties.output          = 'LQ_pars';
IntegralS.calcProperties.model           = 'MKM';
IntegralS.calcProperties.calculusType    = 'rapidMKM_Kase2008_corrected_beta';
%IntegralS.calcProperties.precision       = 0.15;
IntegralS.calcProperties.parallelismType = '0';
IntegralS.calcProperties.cellType        = 'HSG';

IntegralS.calcProperties.modelParam.alpha0   = 0.0;
IntegralS.calcProperties.modelParam.beta0    = 0.05;
IntegralS.calcProperties.modelParam.rNucleus = 3.1;
IntegralS.calcProperties.modelParam.rDomain  = 0.35;

IntegralS.calcProperties.ion = 'H';
IntegralS.calcProperties.trackMode = 'histogram';
IntegralS.calcProperties.energies  = energies;


idx = strfind(IntegralS.wDir, '\');
dirPath = IntegralS.wDir(idx(1)+1:end);
dirPath(strfind(dirPath, '\')) = '/';

IntegralS.wDirWsl = ['/mnt/d/',dirPath];
IntegralS.genParameterFile();
IntegralS.survivalExecutionCommand = ['wsl ',IntegralS.wDirWsl,'/', IntegralS.survivalParameterFileName, '.txt'];
IntegralS.execute();

%% ReadOut
IntegralS.readSingleIonLUT(1);
IntegralS.plotZDSingleIon(0);
