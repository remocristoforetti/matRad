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


pln.machine         = 'newGeneric';
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

%% Override stf
stf = OverrideStf(stf,1);
eIdx = 37;
stf.ray.energy = machine.data(eIdx).energy;
%% DoseCalc

dij_MKM = matRad_calcParticleDose(ct,stf,pln,cst,0);
[MKM_alphaBixel, MKM_betaBixel] = pln.bioParam.calcLQParameter(machine.data(eIdx).depths,machine.data(eIdx),dij_MKM.tissueParameters,[1:size(machine.data(eIdx).depths,1)]);

pln.machine = 'newGenericResampled';
modelName = 'tabulatedRBEModel';
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName,pln.machine);
pln.bioParam.RBEtable = 'RBEtable_rapidMKM_Kase2008_corrected_beta_Survival_rn_0.35';
pln.bioParam.weightMode = 'eD';

dij_MKM_res = matRad_calcParticleDose(ct,stf,pln,cst,0);
[MKM_alphaBixel_res, MKM_betaBixel_res] = pln.bioParam.calcLQParameter(machine.data(eIdx).depths,machine.data(eIdx),dij_MKM_res.tissueParameters,[1:size(machine.data(eIdx).depths,1)]);

%resultGUI_MKM = matRad_calcCubes(1,dij_MKM,1);
% pln.bioParam.RBEtable = 'ExampleMKMInaniwa2010_rn29';
% dij_MKM_rn29 = matRad_calcParticleDose(ct,stf,pln,cst,0);


%% Recalculation LEM
%BaseData LEM
modelName = 'LEM';
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName, pln.machine);
dij_LEM = matRad_calcParticleDose(ct,stf,pln,cst,0);
[LEM_alphaBixel, LEM_betaBixel] = pln.bioParam.calcLQParameter(machine.data(eIdx).depths,machine.data(eIdx),dij_LEM.tissueParameters,[1:size(machine.data(eIdx).depths,1)]);

%resultGUI_LEM = matRAd_calcCubes(1,dij_LEM,1);
%tabulated LEM
modelName = 'tabulatedRBEModel';
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName,pln.machine);
pln.bioParam.weightMode = 'LET';
pln.bioParam.RBEtable = 'RBEtable_rapidLEM_Russo2011_Survival_Dt_60';
dij_LEM_tab = matRad_calcParticleDose(ct,stf,pln,cst,0);
[LEM_alphaBixel_tab, LEM_betaBixel_tab] = pln.bioParam.calcLQParameter(machine.data(eIdx).depths,machine.data(eIdx),dij_LEM_tab.tissueParameters,[1:size(machine.data(eIdx).depths,1)]);
%resultGUI_LEM_tab = matRad_calcCubes(1,dij_LEM_tab,1);

pln.bioParam.weightMode = 'eD';
dij_LEM_tab_eD = matRad_calcParticleDose(ct,stf,pln,cst,0);
[LEM_alphaBixel_tab_eD, LEM_betaBixel_tab] = pln.bioParam.calcLQParameter(machine.data(eIdx).depths,machine.data(eIdx),dij_LEM_tab_eD.tissueParameters,[1:size(machine.data(eIdx).depths,1)]);

%% Plot alpha/beta bixels

figure;
plot(machine.data(eIdx).depths, MKM_alphaBixel, '.-');


hold on;
plot(machine.data(eIdx).depths, LEM_alphaBixel_tab, '.-');
plot(machine.data(eIdx).depths, LEM_alphaBixel, '.-');
plot(machine.data(eIdx).depths, LEM_alphaBixel_tab_eD, '.-');
plot(machine.data(eIdx).depths, MKM_alphaBixel_res, '.-');

grid on;
xlabel('dethps [mm]');
title('alpha bixel');
legend('MKM', 'tabulated LEM', 'base Data LEM', 'tabulated LEM eD', 'MKM resampled');


% figure;
% plot(machine.data(31).depths, MKM_betaBixel, '.-');
% hold on;
% plot(machine.data(31).depths, LEM_betaBixel_tab, '.-');
% plot(machine.data(31).depths, LEM_betaBixel, '.-');
% grid on;
% xlabel('dethps [mm]');
% title('beta bixel');
% legend('MKM', 'tabulated LEM', 'base Data LEM');

%% Multiple energies
stf = OverrideStf(stf,1);
eIdx = [31:39];

MKM_alphaBixel = [];
MKM_betaBixel = [];
modelName = 'tabulatedRBEModel';

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName,pln.machine);

pln.bioParam.RBEtable = 'RBEtable_rapidMKM_Kase2008_corrected_beta_Survival_rn_0.35';
pln.bioParam.weightMode = 'eD';

for k=1:size(eIdx,2)
    stf.ray.energy = machine.data(eIdx(k)).energy;
    dij_MKM = matRad_calcParticleDose(ct,stf,pln,cst,0);
    [alphaBixel, betaBixel] = pln.bioParam.calcLQParameter(machine.data(eIdx(k)).depths,machine.data(eIdx(k)),dij_MKM.tissueParameters,[1:size(machine.data(eIdx(k)).depths,1)]);

    MKM_alphaBixel = [MKM_alphaBixel, {alphaBixel}];
    MKM_betaBixel = [MKM_betaBixel, {betaBixel}];

end

modelName = 'LEM';
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName, pln.machine);

LEM_alphaBixel = [];
LEM_betaBixel  = [];

for k=1:size(eIdx,2)
    stf.ray.energy = machine.data(eIdx).energy;
    dij_LEM = matRad_calcParticleDose(ct,stf,pln,cst,0);
    [alphaBixel, betaBixel] = pln.bioParam.calcLQParameter(machine.data(eIdx(k)).depths,machine.data(eIdx(k)),dij_LEM.tissueParameters,[1:size(machine.data(eIdx(k)).depths,1)]);
    LEM_alphaBixel = [LEM_alphaBixel, {alphaBixel}];
    LEM_betaBixel = [LEM_betaBixel, {betaBixel}];

end

%% Plot multiple alpha bixels
figure;
for k=1:size(eIdx,2)
    plot(machine.data(eIdx(k)).depths, MKM_alphaBixel{k}, '.-');
    hold on;
    plot(machine.data(eIdx(k)).depths, LEM_alphaBixel{k}, '--');
end

grid on;

grid minor;
%% Plot RBExD
MKM_p = matRad_calcCubes(1,dij_MKM,1);
LEM_t_p = matRad_calcCubes(1,dij_LEM_tab,1);
LEM_p = matRad_calcCubes(1,dij_LEM,1);

x = [1:160]*ct.resolution.y - ct.resolution.y/2;

figure;

plot(x, squeeze(sum(MKM_p.RBExD(:,:,:), [2 3])), '.-');
hold on;
plot(x, squeeze(sum(LEM_t_p.RBExD(:,:,:), [2 3])), '.-');
plot(x, squeeze(sum(LEM_p.RBExD(:,:,:), [2 3])), '.-');
grid on;
xlabel('depth [mm]');
legend('MKM', 'LEM tabulated', 'LEM base Data');

% figure;
% %plot(machine.data(eIdx).depths, MKM_alphaBixel, '.-');
% plot(x, squeeze(sum(MKM_p.RBExD(:,:,:), [2 3])), '.-');
% hold on;
% plot(x, squeeze(sum(MKM_p.physicalDose(:,:,:), [2 3])), '.-');
% legend('RBExD', 'physical Dose');
% grid on;

%% Plot RBEtables and alpha/LET

load('RBEtable_rapidLEM_Russo2011_Survival_Dt_60.mat');
RBEtable_LEM = RBEtable;

load('RBEtable_rapidMKM_Kase2008_corrected_beta_Survival_rn_0.35.mat');
RBEtable_MKM = RBEtable;

eneMKM = RBEtable_MKM.data(1).energies;
eneLEM = RBEtable_LEM.data(1).energies;


ion = 6;

alphaMKM = RBEtable_MKM.data(1).alpha(:,ion);
alphaLEM = RBEtable_LEM.data(1).alpha(:,ion);
letsMKM = arrayfun(@(e) pln.bioParam.computeLET(e,ion), eneMKM);
letsLEM = arrayfun(@(e) pln.bioParam.computeLET(e,ion), eneLEM);
betaMKM = RBEtable_MKM.data(1).beta(:,ion);
betaLEM = RBEtable_LEM.data(1).beta(:,ion);

figure;
semilogx(eneMKM, alphaMKM, '.-');

hold on;
semilogx(eneLEM, alphaLEM, '.-');
grid on;
xlabel('E [MeV/u]');
ylabel('alpha [Gy^{-1}]');
legend('MKM', 'LEM');

% figure;
% semilogx(letsMKM, alphaMKM, '.-');
% hold on;
% semilogx(letsLEM, alphaLEM, '.-');
% grid on;
% xlabel('LET [keV/um]');
% ylabel('alpha [Gy^{-1}]');
% legend('MKM', 'LEM')

figure;
semilogx(eneMKM, betaMKM, '.-');
hold on;
semilogx(eneLEM, betaLEM, '.-');
grid on;
xlabel('E [MeV/u]');
ylabel('beta [Gy^{-2}]');
legend('MKM', 'LEM');
% 
% figure;
% semilogx(letsMKM, betaMKM, '.-');
% hold on;
% semilogx(letsLEM, betaLEM, '.-');
% grid on;
% xlabel('LET [keV/um]');
% ylabel('beta [Gy^{-1}]');
% legend('MKM', 'LEM')


%% alphaBixel betaBixel, this is exactly the same as computing the varAlpha/varBeta

filenames = dir([matRad_cfg.matRadRoot, filesep, 'bioModels', filesep, 'RBEtables', filesep, 'RBEtable_rapidMKM_Kase2008_corrected_beta_Survival_rn_*.mat']);
rn = [];
for k=1:size(filenames,1)

    modelName = 'tabulatedRBEModel';
    pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName,pln.machine);
    
    idxDot = find(filenames(k).name == '.',1,'last');
    nameTable = filenames(k).name(1:idxDot-1);
    pln.bioParam.RBEtable = nameTable;

    dij_MKM = matRad_calcParticleDose(ct,stf,pln,cst,0);
    
    %MKM_p(:,k) = sum(reshape(full(dij_MKM.mAlphaDose{1}), [dij_MKM.doseGrid.dimensions]), [2 3]);

    [MKM_alhpaBixel(:,k), MKM_betaBixel(:,k)] = pln.bioParam.calcLQParameter(machine.data(31).depths,machine.data(31),dij_MKM.TissueParameters,[1:size(machine.data(31).depths,1)]);
    
    load([filenames(k).folder, filesep, filenames(k).name]);
    rn = [rn, RBEtable.meta.modelParameters.rDomain];
end
    
%     pln.bioParam.RBEtable = 'ExampleMKMInaniwa2010_rn29';
% [MKM_alhpaBixel_rn29, MKM_betaBixel_rn29] = pln.bioParam.calcLQParameter(machine.data(31).depths,machine.data(31),dij_MKM_rn29.TissueParameters,[1:size(machine.data(31).depths,1)]);


modelName = 'tabulatedRBEModel';
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName,pln.machine);
pln.bioParam.weightMode = 'LET';
pln.bioParam.RBEtable = 'RBEtable_rapidLEM_Russo2011_Survival_Dt_30';

dij_LEM_tab = matRad_calcParticleDose(ct,stf,pln,cst,0);
[LEM_alphaBixel_tab, LEM_betaBixel_tab] = pln.bioParam.calcLQParameter(machine.data(31).depths,machine.data(31),dij_LEM_tab.TissueParameters,[1:size(machine.data(31).depths,1)]);

modelName = 'LEM';
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName, pln.machine);
[LEM_alhpaBixel, LEM_betaBixel] = pln.bioParam.calcLQParameter(machine.data(31).depths,machine.data(31),dij_LEM.TissueParameters,[1:size(machine.data(31).depths,1)]);

%% Plot bixelAplha bixelbeta
figure;

lege = [];
for k=6%:size(filenames,1)

    plot(machine.data(31).depths, MKM_alhpaBixel(:,k), '.-');
    hold on;
    lege = [lege, {['MKM rn = ', num2str(rn(k))]}];
end
%plot(machine.data(31).depths, MKM_alhpaBixel_rn29, '.-');
plot(machine.data(31).depths,LEM_alphaBixel_tab,'.-', 'color', 'r');


plot(machine.data(31).depths,LEM_alhpaBixel, '.-', 'color', 'k');
grid on;

legend([lege,'LEM tabulated', 'LEM']);
title('alphaBixel');

%MKM_p = sum(reshape(full(dij_MKM.mAlphaDose{1}), [dij_MKM.doseGrid.dimensions]), [2 3]);
%MKM_rn29_p = sum(reshape(full(dij_MKM_rn29.mAlphaDose{1}), [dij_MKM_rn29.doseGrid.dimensions]), [2 3]);
%LEM_p = sum(reshape(full(dij_LEM.mAlphaDose{1}), [dij_LEM.doseGrid.dimensions]), [2 3]);

% figure;
% plot([1:96],MKM_p, '.-');
% hold on;
% plot([1:96],MKM_rn29_p, '.-');
% plot([1:96],LEM_p, '.-');
% grid on;
% legend('MKM rn = 0.35', 'MKM rn = 0.29', 'LEM');

%% Plot varAlpha and beta

MKM_alpha = dij_MKM.TissueParameters.baseParam.varAlpha;
MKM_beta =  dij_MKM.TissueParameters.baseParam.varBeta;

LEM_alpha = machine.data(31).alpha;
LEM_beta = machine.data(31).beta;

figure;
plot(dij_MKM.TissueParameters.baseParam.d, MKM_alpha(:,1), '.-');
hold on;
plot(machine.data(31).depths, LEM_alpha(:,1), '.-');
grid on;
xlabel('depth [mm]');
ylabel('alpha [Gy^{-1}]');
legend('MKM', 'LEM')

figure;
plot(dij_MKM.TissueParameters.baseParam.d, MKM_beta(:,1), '.-');
hold on;

plot(machine.data(31).depths, LEM_beta(:,1), '.-');
grid on;
xlabel('depth [mm]');
ylabel('alpha [Gy^{-1}]');
legend('MKM', 'LEM');

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


energies = logspace(-1,log10(300),500);
IntegralS = integralDoseSurvival();
IntegralS.wDir = [matRad_cfg.matRadRoot, filesep, 'survival'];

%This is source path see from wsl
IntegralS.survivalSourcePath = '/mnt/c/Users/r408i/Desktop/r408i_data/Survival';

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

%this is working directory from the wsl perspective
IntegralS.wDirWsl = ['/mnt/c/',dirPath];
IntegralS.calcProperties.trackMode = 'histogram';
IntegralS.calcProperties.energies  = energies;

idx = strfind(IntegralS.wDir, '\');
dirPath = IntegralS.wDir(idx(1)+1:end);
dirPath(strfind(dirPath, '\')) = '/';


ions  = {'H', 'He', 'Li', 'Be', 'B', 'C'};
ionsMassNumber = [1 4 7 9 11 12]; 
for k=1:size(ions,2)
    IntegralS.calcProperties.energies  = energies.*ionsMassNumber(k);
    IntegralS.calcProperties.ion = ions{k};
    IntegralS.genParameterFile();
    IntegralS.survivalExecutionCommand = ['wsl ',IntegralS.wDirWsl,'/', IntegralS.survivalParameterFileName, '.txt'];

    %Execution with system does not work, check -> Solved, if does not work
    %again, check that default destribution is Ubuntu and not Ubunutu-20.04.
    %Else, from prompt wsl --setdefault Ubuntu
    matRad_cfg.dispInfo('Executing Survival...');
    a = IntegralS.execute();
    if a == 0
        matRad_cfg.dispInfo('done. \n');
    else
        matRad_cfg.dispInfo('error. \n');
    end
end
%% ReadOut

%these alphaE and betaE are E,I specific. in case of MKM, this is zD(E,I),
%need to do (alphaX + betaX*zD(E,I)) and then couple to fluence, in order
%to get varAlpha varBeta
[alphaE, betaE] = IntegralS.readMultipleIonLUT(ions);
%IntegralS.plotZDSingleIon(0);
figure;

for k=1:size(ions,2)
    semilogx(energies, alphaE{k}, '.-');
    hold on;
end