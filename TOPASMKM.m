matRad_rc;
matRad_cfg = MatRad_Config.instance();
load('BOXPHANTOM.mat');

%% Setup

pln.numOfFractions        = 25;
pln.propStf.gantryAngles  = 0;
pln.propStf.couchAngles   = 0;
pln.propStf.bixelWidth    = 3;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;


pln.radiationMode   = 'protons';

pln.machine         = 'Generic';
pln.multScen = matRad_multScen(ct,'nomScen');

pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]
stf = matRad_generateStf(ct,cst,pln);

%stf = OverrideStf(stf);

load([pln.radiationMode, '_',pln.machine]);
quantityOpt = 'RBExD';

modelName = 'MKM';
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);

%% Dose calculation

% load('D:\matRad_gitHubRemo_MKM_TOPAS\MKM\Spectrum_ProtonsMultiEnergy_Beta005.mat');
% [~,Ex] = intersect([machine.data(:).energy],stf.ray(1).energy);
% for k=1:size(Ex,1)
%      
%        machine.data(Ex(k)).Spe = Spe(k);
% end
   
if (strcmp(modelName, 'MKM') && ~strcmp(pln.machine,'GenericMKM'))
   try
      load([pln.radiationMode, '_', pln.machine, 'MKM_Beta005']);
      pln.machine = 'GenericMKM_Beta005';
   catch
      matRad_cfg.dispWarning('Unable to find machine data');
   end
end

AvailableDataIdx = find(arrayfun(@(x) ~isempty(x.Spe), [machine.data]));
if ~size(intersect([machine.data(AvailableDataIdx).energy],unique([stf.ray(:).energy])),2) == size(unique(stf.ray(1).energy),2)
   matRad_cfg.dispWarning('Energies are not correct');
end

dij_MKM = matRad_calcParticleDose(ct,stf,pln,cst,0);

resultGUI = matRad_fluenceOptimization(dij_MKM,cst,pln);

%% Recompute results
quantityOpt = 'RBExD';
modelName = 'MCN';

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);
dij_MCN = matRad_calcParticleDose(ct,stf,pln,cst,0);

RBExD_MCN = matRad_calcCubes(resultGUI.w,dij_MCN,1);

quantityOpt = 'RBExD';
modelName = 'WED';
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);
dij_WED = matRad_calcParticleDose(ct,stf,pln,cst,0);

RBExD_WED = matRad_calcCubes(resultGUI.w,dij_WED,1);

%% Plot Results
Profile_MKM = resultGUI.RBExD(:,80,80);
Profile_MCN = RBExD_MCN.RBExD(:,80,80);
Profile_WED = RBExD_WED.RBExD(:,80,80);

x = [1:160]*3 - 3/2;

figure;
plot(x, Profile_MCN, '.-');
hold on;
plot(x, Profile_WED, '.-');
plot(x, Profile_MKM, '.-');
grid on;
legend('MCN','WED', 'MKM');
%% Topas Setup

stf = OverrideStf(stf);
topas_cfg = MatRad_TopasConfig;
topas_cfg.numOfRuns = 1;
numHistoriesPerBixel = 100000;
topas_cfg.numHistories = ceil(numHistoriesPerBixel*size(stf.ray.energy,2));
topas_cfg.outputType = 'DICOM';
topas_cfg.scoreDose = false;
topas_cfg.addVolumeScorers = true;


topasBaseData = MatRad_TopasBaseData(machine,stf);
w = (1/size(stf.ray.energy,2))*ones(size(stf.ray.energy,2),1);ls

% ctT = OverrideCT(ct);
% stf.isoCenter = [ctT.resolution.x/2 ctT.resolution.y/2, 0.5];
ctT = ct;

ctT.cubeHU = {zeros(ctT.cubeDim)};

topas_cfg.writeAllFiles(ctT,pln,stf,topasBaseData,w);
WriteSplitterFeature(stf, topas_cfg.workingDir);

pathLUT = 'X:\MKMTEST\Proton.mat';
WriteProtonLUT(pathLUT, topas_cfg.workingDir);
%% ReadOut TopasData
filenames = dir('D:\matRad_gitHubRemo_MKM_TOPAS\MKM\MKM_SimulationData_Beta005\*.dcm');
zMix = [];
d = [1:480] - 0.5; %[1:160]*3 -3/2;%-0.5;
for k=1:size(filenames)
    infos = dicominfo(strcat(filenames(1).folder, filesep, filenames(k).name));
    zMix = [zMix, double(dicomread(strcat(filenames(1).folder, filesep, filenames(k).name)))*infos.DoseGridScaling];
end

filenames = dir('X:\MKMTEST\Results\LUTProton_Seed_1\*.dcm');
zMix_1 = [];
for k=1:size(filenames)
    zMix_1 = [zMix_1, dicomread(strcat(filenames(1).folder, filesep, filenames(k).name))];

end

%toplot = [1:7];
%figure;
% plot(d, zMix_1(:,toplot), '.-');
% hold on;
%plot(d, zMix(:,toplot), '.-');


for k = 1:size(zMix,2)
    Spe(k).d = d;
    Spe(k).zMix = zMix(:,k);
end
%% Comparison of REBxD
quantityOpt = 'RBExD';
modelName = 'MCN';
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);
dij_MCN = matRad_calcParticleDose(ct,stf,pln,cst,0);

quantityOpt = 'RBExD';
modelName = 'WED';
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);
dij_WED = matRad_calcParticleDose(ct,stf,pln,cst,0);

modelName = 'MKM';
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);
dij_MKM = matRad_calcParticleDose(ct,stf,pln,cst,0);
%% Cubes

%weights = ones(dij_MCN.totalNumOfBixels,1);
weights = exp([1:23])'.*ones(dij_MCN.totalNumOfBixels,1)/(3*10^7);

e_MCN = dij_MCN.mAlphaDose{1}*weights + (dij_MCN.mSqrtBetaDose{1}*weights).^2;
EffectCube_MCN = reshape(full(e_MCN), dij_MCN.doseGrid.dimensions);
Cube_pD = reshape(full(dij_MCN.physicalDose{1}*weights), dij_MCN.doseGrid.dimensions);

e_WED = dij_WED.mAlphaDose{1}*weights + (dij_WED.mSqrtBetaDose{1}*weights).^2;
EffectCube_WED = reshape(full(e_WED), dij_WED.doseGrid.dimensions);

e_MKM = dij_MKM.mAlphaDose{1}*weights + (dij_MKM.mSqrtBetaDose{1}*weights).^2;
EffectCube_MKM = reshape(full(e_MKM), dij_MKM.doseGrid.dimensions);


%% Plots
Central = 239;
Profile_MCN_effect = EffectCube_MCN(:,Central,Central);
Profile_WED_effect = EffectCube_WED(:,Central,Central);
Profile_MKM_effect = EffectCube_MKM(:,Central,Central);
Profile_pD = Cube_pD(:,Central,Central);

AlphaX = cst{1,5}.alphaX;
BetaX = cst{1,5}.betaX;

RBExD_Cube_MCN = sqrt(EffectCube_MCN/BetaX + (AlphaX/(2*BetaX))^2) - AlphaX/(2*BetaX);
RBExD_Cube_WED = sqrt(EffectCube_WED/BetaX + (AlphaX/(2*BetaX))^2) - AlphaX/(2*BetaX);
RBExD_Cube_MKM = sqrt(EffectCube_MKM/BetaX + (AlphaX/(2*BetaX))^2) - AlphaX/(2*BetaX);

Profile_MCN_RBExD = RBExD_Cube_MCN(:,Central,Central);
Profile_WED_RBExD = RBExD_Cube_WED(:,Central,Central);
Profile_MKM_RBExD = RBExD_Cube_MKM(:,Central,Central);

d = [1:478]*pln.propDoseCalc.doseGrid.resolution.x;

figure;
plot(d,Profile_MCN_RBExD,'.-');
hold on;
plot(d,Profile_WED_RBExD,'.-');
hold on;
plot(d, Profile_MKM_RBExD, '.-');
plot(d, Profile_pD, '.-', 'color', 'k');
legend('MCN', 'WED','MKM', 'physical Dose');

grid on;
%xlim([100 350]);
% figure;
% imagesc((EffectCube_MCN(:,:,80) - EffectCube_MKM(:,:,80))./(EffectCube_MCN(:,:,80)));
% colorbar;

%% ReadOut Spectra
% DataFiles = dir(fullfile(regexprep(topas_cfg.workingDir,'\MCrun','\Results'), '*.csv'));
% BinLimits = [0 150];
% nBins = 150;
% BinWidth = (BinLimits(2)-BinLimits(1))/nBins;
% E = linspace(BinLimits(1) + BinWidth/2, BinLimits(2)-BinWidth/2,nBins);
% 
% for k=1:1%size(DataFiles,1)
%     %Data = import_csv_Results(strcat(DataFiles(k).folder, '\',DataFiles(k).name),[11 Inf],nBins+6);
%     Data = csvread(strcat(DataFiles(k).folder, '\',DataFiles(k).name),10, 3);%, [10,3,109, nBins+3]);
%     Spectrum(k).depths = Data(2:end,1); %Data(:,1);
%     Spectrum(k).Energies = Data(2:end,1:end); %Data(:,2:end);
% end
% 
% 
% figure;
% plot(E,Spectrum(1).Energies(1,:),'.-');
% % figure;
% % histogram(Spectrum(1).Energies(end,:));
% hold on;
% plot(E,Spectrum(1).Energies(30,:),'.-');
% % hold on;