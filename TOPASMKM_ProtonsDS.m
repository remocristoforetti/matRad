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


pln.radiationMode   = 'protons';

pln.machine         = 'generic_TOPAS_mcSquare';
pln.multScen = matRad_multScen(ct,'nomScen');

pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]
stf = matRad_generateStf(ct,cst,pln);
pln.propDoseCalc.calcLET = 1;
%stf = OverrideStf(stf);

load([pln.radiationMode, '_',pln.machine]);
quantityOpt = 'RBExD';

modelName = 'MKM';
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);

%% Dose calculation
% 
% load('D:\matRad_gitHubRemo_MKM_TOPAS\MKM\Spectrum_CarbonsMultiEnergy.mat');
% [~,Ex] = intersect([machine.data(:).energy],stf.ray(1).energy);
% for k=1:size(Ex,1)
%      
%        machine.data(Ex(k)).Spe = Spe(k);
% end
   
if (strcmp(modelName, 'MKM') && ~strcmp(pln.machine,'GenericMKM_DS_rn25'))
   try
      load([pln.radiationMode, '_', pln.machine, 'MKM_DS_rn25']);
      pln.machine = 'GenericMKM_DS_rn25';
   catch
      matRad_cfg.dispWarning('Unable to find machine data');
   end
end

AvailableDataIdx = find(arrayfun(@(x) ~isempty(x.Spe), [machine.data]));
if ~size(intersect([machine.data(AvailableDataIdx).energy],unique([stf.ray(:).energy])),2) == size(unique(stf.ray(1).energy),2)
   matRad_cfg.dispWarning('Energies are not correct');
end

pln.machine = 'GenericMKM_DS_rn25';
dij_MKM = matRad_calcParticleDose(ct,stf,pln,cst,0);
resultGUI = matRad_fluenceOptimization(dij_MKM,cst,pln);

%% Recompute results
quantityOpt = 'RBExD';
modelName = 'MCN';
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);

pln.machine = 'Generic';
dij_MCN = matRad_calcParticleDose(ct,stf,pln,cst,0);
RBExD_MCN = matRad_calcCubes(resultGUI.w,dij_MCN,1);

modelName = 'WED';
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);

pln.machine = 'Generic';
dij_WED = matRad_calcParticleDose(ct,stf,pln,cst,0);

RBExD_WED = matRad_calcCubes(resultGUI.w,dij_WED,1);

%% Plot Results
Profile_MKM = resultGUI.RBExD(:,80,80);
Profile_MCN = RBExD_MCN.RBExD(:,80,80);
Profile_WED = RBExD_WED.RBExD(:,80,80);

Profile_pd = resultGUI.physicalDose(:,80,80);
Profile_MKM_LET = resultGUI.LET(:,80,80);

x = [1:160]*3 - 3/2;

figure;
subplot(3,1,[1,2])
plot(x, Profile_MCN, '.-');
hold on;
plot(x, Profile_MKM, '.-');
plot(x, Profile_WED, '.-');
plot(x, Profile_pd, '.-');
ylabel('Dose [Gy]');
yyaxis right;
plot(x, Profile_MKM_LET, '.-');
grid on;
ylabel('KeV/\mu');
xlabel('Depth [mm]');
legend('MCN', 'MKM', 'WED', 'physical Dose', 'LET');
subplot(3,1,3)
plot(x, (Profile_MKM - Profile_MCN)./Profile_MKM, '.-');
hold on;
plot(x, (Profile_MKM - Profile_WED)./Profile_MKM, '.-', 'color', [0.9290, 0.6940, 0.1250]);
legend('(MKM - MCN)/MKM', '(MKM - WED)/MKM', 'Location','SouthWest');
grid on;
xlim([100 280]);
%% Topas Setup

stf = OverrideStf(stf);

topas_cfg = matRad_TopasConfig;

topas_cfg.scorer.filename = 'TOPAS_scorer_ZMixMultiIon.txt.in';
topas_cfg.numOfRuns = 1;
numHistoriesPerBixel = 100000;
topas_cfg.numHistories = ceil(numHistoriesPerBixel*size(stf.ray.energy,2));

topas_cfg.scorer.outputType = 'DICOM';
topas_cfg.scorer.reportQuantity = {'Sum'};

topas_cfg.scorer.doseToMedium = 0;
w = (1/size(stf.ray.energy,2))*ones(size(stf.ray.energy,2),1);

% ctT = OverrideCT(ct);
% stf.isoCenter = [ctT.resolution.x/2 ctT.resolution.y/2, 0.5];
ctT = ct;


ctT.cubeHU = {zeros(ctT.cubeDim)};
ctT.ctGrid.resolution = ct.resolution;

%% Read data example
Data = dicomread('D:\matRad_gitHubRemo_MKM_TOPAS\Ds.dcm');
Data = squeeze(Data);
DataP = permute(Data, [2,3,1]);
depthValues = [1:160]*3 -3/2;
EMin = 0.1;
EMax = 130;
nEBins = 1300;

EBinWidth = (EMax-EMin)/nEBins;

EdMin = 0;
EdMax = 3;
nEdBins = 300;

EdBinWidth = (EdMax-EdMin)/nEdBins;

EBins = linspace(EMin,EMax,nEBins)+EBinWidth/2;
EdepBins = linspace(EdMin,EdMax,nEdBins)+ EdBinWidth/2;

doubleS = DS_protons();
doubleS.d = depthValues;
doubleS.E = EBins;
doubleS.Edep = EdepBins;
doubleS.rn = 0.25;
doubleS.RN = 3.9;
doubleS.Beta_Tissue = 0.05;

RawD = DataP(2:end,2:end,:);
doubleS.Phi = double(RawD);

doubleS.computeZD();

% Profile_Zmix = dicomread('D:\matRad_gitHubRemo_MKM_TOPAS\ZMixOutput.dcm');
% infos = dicominfo('D:\matRad_gitHubRemo_MKM_TOPAS\ZMixOutput.dcm');
% 
% Profile_Zmix = double(Profile_Zmix)*infos.DoseGridScaling;

%figure;
%plot(doubleS.d, doubleS.zD, '.-');

%plot(doubleS.d, Profile_Zmix, '.-');

%xlim([0,110]);


%% Read Double Spectrum example
v = dicomread('D:\matRad_gitHubRemo_MKM_TOPAS\DS_proton_ResultsMultienergy\Ion_1\Ds_Run_0000.dcm');%dicomread('D:\matRad_gitHubRemo_MKM_TOPAS\DS_Fragments\ZMixOutput_C.dcm');%dicomread('D:\matRad_gitHubRemo_MKM_TOPAS\ZMixOutput_C.dcm');%dicomread('D:\matRad_gitHubRemo_MKM_TOPAS\Ds.dcm');
%infos = dicominfo('D:\matRad_gitHubRemo_MKM_TOPAS\Ds.dcm');
v = double(squeeze(v));%*infos.DoseGridScaling;

m = permute(v, [2,3,1]);
EMax = 130;
EMin = 0.1;
EBins = 1300;

EdMax = 3;
EdMin = 0;
EdBins = 300;


%% Plot example
toplot = [1:10:110];

for k=toplot
   figure;
   imagesc(double(m(2:end,2:end,k))./sum(m, [1,2,3]));
   title('Z = ',num2str(k - 0.5));
   currXtick = get(gca, 'XTick');
   currYtick = get(gca, 'YTick');

   
   set(gca,'XTick',currXtick, 'XTickLabel', linspace(currXtick(1)*EMax/(currXtick(end)+currXtick(1)/2),EMax,size(currXtick,2)));

   set(gca,'YTick', currYtick, 'YTickLabel', linspace(currYtick(1)*EdMax/(currYtick(end) + currYtick(1)/2),EdMax,size(currYtick,2)));
   xlabel('E [MeV]');
   ylabel('Edep [MeV]');

   colorbar;
end

%% Read Multienergy
doubleS = DS_MultiIon();

doubleS.wDir = 'D:\matRad_gitHubRemo_MKM_TOPAS\DS_proton_ResultsMultienergy';
doubleS.d = [1:480] - 0.5;%[1:240]*2 - 1;
doubleS.LUT = 'D:\matRad_gitHubRemo_MKM_TOPAS\MKM\LUTs_rn25_RN_39_Beta005';
doubleS.rn = 0.25;
doubleS.RN = 3.9;
doubleS.Beta_Tissue = 0.05;

doubleS.Ions = {'protons'};

for k=1:size(doubleS.Ions,2)
   doubleS.EParam(k).EMax = 130.1;
   doubleS.EParam(k).EMin = 0.1;
   doubleS.EParam(k).nEBins = 1300;
   
   doubleS.EParam(k).EdMax = 3;
   doubleS.EParam(k).EdMin = 0;
   doubleS.EParam(k).nEdBins = 300;
end

doubleS.computeBinning();

for k=1:23
   doubleS.analyzeRawData(k-1);
   doubleS.computeZD();

   Spe(k).d = doubleS.d;
   Spe(k).zMix  = doubleS.zD;
   display(['Computing zD for Energy: ', num2str(k+1)]);

end
%
figure;
for k=1:10
   plot(Spe(k).d, Spe(k).zMix, '.-');
   hold on;
end


% for k=1:23
%    Spe(k).d = tmp(k).d;
%    Spe(k).zMix = tmp(k).zD;
% end
%% Generate machine data
[~,Ex] = intersect([machine.data(:).energy],stf.ray(1).energy);
for k=1:size(Ex,1)
     
       machine.data(Ex(k)).Spe = Spe(k);
end

%% profile comparison for basedata different LUTs
OldBaseData = load('protons_GenericMKM');
OldProfile = OldBaseData.machine.data(28).Spe;

figure;
plot(Spe(1).d, Spe(1).zD, '.-');
hold on;
plot(OldProfile.d, OldProfile.zMix, '.-');
grid on;
legend('New LUT', 'Old LUT');

%% Profile comparison MKM vs MCN vs WED
M = load(['protons_GenericMKM_DS_rn25.mat']);
M = M.machine;


EIdx = 28+22;
depths = M.data(EIdx).depths;

MCN_Param = matRad_bioModel('protons', 'RBExD', 'MCN');
MCN_Profile = MCN_Param.calcLQParameter(depths,M.data(EIdx),ones(size(depths)),0.1*ones(size(depths)),0.05*ones(size(depths)),2*ones(size(depths)));

WED_Param = matRad_bioModel('protons', 'RBExD', 'WED');
WED_Profile = WED_Param.calcLQParameter(depths,M.data(EIdx),ones(size(depths)),0.1*ones(size(depths)),0.05*ones(size(depths)),2*ones(size(depths)));


MKM_Param = matRad_bioModel('protons', 'RBExD', 'MKM');
MKM_Profile = MKM_Param.calcLQParameter(depths,M.data(EIdx),ones(size(depths)),0.1*ones(size(depths)),0.05*ones(size(depths)),2*ones(size(depths)));


figure;
plot(depths, MCN_Profile, '.-');
hold on;
plot(depths, MKM_Profile, '.-');
plot(depths, WED_Profile, '.-');

grid on;
legend('MCN', 'MKM', 'WED');