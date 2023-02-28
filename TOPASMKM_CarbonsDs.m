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
%stf = OverrideStf(stf);
load([pln.radiationMode, '_',pln.machine]);
quantityOpt = 'RBExD';

modelName = 'MKM';
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);
% stf = OverrideStf(stf);
% 
% 
% dij = matRad_calcParticleDose(ct, stf, pln, cst, 0);
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

%pln.machine = 'GenericMKM';
pln.propDoseCalc.airOffsetCorrection = true;
dij_MKM = matRad_calcParticleDose(ct,stf,pln,cst,0);
resultGUI = matRad_fluenceOptimization(dij_MKM,cst,pln);

%% Recompute results

quantityOpt = 'RBExD';
modelName = 'LEM';
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);

pln.machine = 'Generic';
dij_LEM = matRad_calcParticleDose(ct,stf,pln,cst,0);

RBExD_LEM = matRad_calcCubes(resultGUI.w,dij_LEM,1);
%% Plot Results
Profile_MKM = resultGUI.RBExD(:,80,80);
Profile_LEM = RBExD_LEM.RBExD(:,80,80);
Profile_pd = resultGUI.physicalDose(:,80,80);
Profile_MKM_LET = resultGUI.LET(:,80,80);

x = [1:160]*3 - 3/2;

figure;
subplot(3,1,[1,2])
plot(x, Profile_LEM, '.-');
hold on;
plot(x, Profile_MKM, '.-');
plot(x, Profile_pd, '.-');
ylabel('Dose [Gy]');
yyaxis right;
%plot(x, Profile_MKM_LET, '.-');
grid on;
ylabel('KeV/\mu');
xlabel('Depth [mm]');
legend('LEM', 'MKM', 'physical Dose');%, 'LET');
subplot(3,1,3)
plot(x, (Profile_MKM - Profile_LEM)./Profile_MKM, '.-');
grid on;
legend('(MKM -LEM)/MKM', 'Location', 'SouthWest');
%% Topas Setup

stf = OverrideStf(stf);

topas_cfg = matRad_TopasConfig;

topas_cfg.numOfRuns = 1;
numHistoriesPerBixel = 100000;
topas_cfg.numHistories = ceil(numHistoriesPerBixel*size(stf.ray.energy,2));

topas_cfg.scorer.outputType = 'DICOM';
topas_cfg.scorer.reportQuantity = {'Sum'};
%topas_cfg.scorer.filename = 'TOPAS_scorer_ZMixMultiIon.txt.in';

%topas_cfg.useOrigBaseData = true;

topas_cfg.scorer.doseToMedium = 0;
w = (1/size(stf.ray.energy,2))*ones(size(stf.ray.energy,2),1);

% ctT = OverrideCT(ct);
% stf.isoCenter = [ctT.resolution.x/2 ctT.resolution.y/2, 0.5];
ctT = ct;

ctT.cubeHU = {zeros(ctT.cubeDim)};

ctT.ctGrid.resolution = ct.resolution;


topas_cfg.writeAllFiles(ct,cst,pln,stf,machine,w);

%% Read data
doubleS = DS_MultiIon();

doubleS.wDir = 'D:\matRad_gitHubRemo_MKM_TOPAS\DS_carbon_MultiEnergyLocal\LUTResult\';
doubleS.d = [1:160]*3 - 3/2;%[1:240]*2 - 1;

doubleS.rn = 0.25;
doubleS.RN = 3.9;
doubleS.Beta_Tissue = 0.05;


%Protons
doubleS.EParam(1).EMax = 250.1;
doubleS.EParam(1).EMin = 0.1;
doubleS.EParam(1).nEBins = 500;

doubleS.EParam(1).EdMax = 12;
doubleS.EParam(1).EdMin = 0;
doubleS.EParam(1).nEdBins = 240;

%He
doubleS.EParam(2).EMax = 250.1;
doubleS.EParam(2).EMin = 0.1;
doubleS.EParam(2).nEBins = 500;

doubleS.EParam(2).EdMax = 12;
doubleS.EParam(2).EdMin = 0;
doubleS.EParam(2).nEdBins = 240;

%Li
doubleS.EParam(3).EMax = 250.1;
doubleS.EParam(3).EMin = 0.1;
doubleS.EParam(3).nEBins = 500;

doubleS.EParam(3).EdMax = 12;
doubleS.EParam(3).EdMin = 0;
doubleS.EParam(3).nEdBins = 240;

%Be
doubleS.EParam(4).EMax = 250.1;
doubleS.EParam(4).EMin = 0.1;
doubleS.EParam(4).nEBins = 500;

doubleS.EParam(4).EdMax = 12;
doubleS.EParam(4).EdMin = 0;
doubleS.EParam(4).nEdBins = 240;

%B
doubleS.EParam(5).EMax = 250.1;
doubleS.EParam(5).EMin = 0.1;
doubleS.EParam(5).nEBins = 500;

doubleS.EParam(5).EdMax = 12;
doubleS.EParam(5).EdMin = 0;
doubleS.EParam(5).nEdBins = 240;

%C
doubleS.EParam(6).EMax = 250.1;
doubleS.EParam(6).EMin = 0.1;
doubleS.EParam(6).nEBins = 500;

doubleS.EParam(6).EdMax = 12;
doubleS.EParam(6).EdMin = 0;
doubleS.EParam(6).nEdBins = 240;

doubleS.computeBinning();

doubleS.Ions = {'protons', 'He', 'Li', 'Be', 'B', 'C'};

%Energies = [2460.031242 2489.177179 2523.404295 2555.041491 2584.251617 2645.318850 2677.054488 2705.875518 2734.448317 2769.476893 2796.254945 2829.173914 2886.299737 2912.624071 2944.570220 2972.973280 2998.839371 3026.140081 3058.427119 3086.921845 3140.908146 3166.921755 3196.787759 3223.982853 3278.619636 3306.112261 3358.625818 3406.661471];

for k=1:size(doubleS.EParam,2)
   doubleS.EParam(k).EMax = 250.1;
   doubleS.EParam(k).EMin = 0.1;
   doubleS.EParam(k).nEBins = 250;
   
   doubleS.EParam(k).EdMax = 12;
   doubleS.EParam(k).EdMin = 0;
   doubleS.EParam(k).nEdBis = 240;
end

doubleS.computeBinning();

for k=1:size(Energies,2)
   doubleS.analyzeRawData(k-1);
   doubleS.computeZD();

   Spe(k).E = Energies(k);
   Spe(k).d = doubleS.d;
   Spe(k).zD  = doubleS.zD;
   display(['Computing zD for Energy: ', num2str(k+1)]);

end

figure;
for k=1:10
   plot(Spe(k).d, Spe(k).zD, '.-');
   hold on;
end
grid on;

% doubleS.computeZD();
% doubleS.computeZD_SingleIons();
%% Plotting
figure;

plot(doubleS.d, doubleS.zD, '.-');
hold on;
for k=1:size(doubleS.Ions,2)
   plot(doubleS.d, doubleS.zD_SingleIons{k}, '.-');
   hold on;
end
legend([{'Total zD'}, Ions]);
grid on;
xlim([0, 130]);
%% Read Double Spectrum
v = dicomread('D:\matRad_gitHubRemo_MKM_TOPAS\DS_carbon_MultiEnergyLocal\LUTResult\Ion_6\ZMixOutput_C_Run_0000.dcm');%dicomread('D:\matRad_gitHubRemo_MKM_TOPAS\DS_Fragments\ZMixOutput_C.dcm');%dicomread('D:\matRad_gitHubRemo_MKM_TOPAS\ZMixOutput_C.dcm');%dicomread('D:\matRad_gitHubRemo_MKM_TOPAS\Ds.dcm');
%infos = dicominfo('D:\matRad_gitHubRemo_MKM_TOPAS\Ds.dcm');
v = double(squeeze(v));%*infos.DoseGridScaling;
m= permute(v, [2,3,1]);

EMax = 500/2;
EMin = 0.1;
EBins = 500;

EdMax = 12;
EdMin = 0;
EdBins = 240;

toplot = [1:5:35];
for k=toplot
   figure;
   imagesc(double(m(2:end,2:end,k))./sum(m, [1,2,3]));
   title('Z = ',num2str(k*2 - 1));
   currXtick = get(gca, 'XTick');
   currYtick = get(gca, 'YTick');

   
   set(gca,'XTick',currXtick, 'XTickLabel', linspace(currXtick(1)*EMax/(currXtick(end)+currXtick(1)/2),EMax,size(currXtick,2)));

   set(gca,'YTick', currYtick, 'YTickLabel', linspace(currYtick(1)*EdMax/(currYtick(end)+currYtick(1)/2),EdMax,size(currYtick,2)));
   xlabel('E [MeV]');
   ylabel('Edep [MeV]');

   colorbar;
   
end
%% Double Spectrum LUTs comparison

doubleS = DS_MultiIon();

doubleS.wDir = 'D:\matRad_gitHubRemo_MKM_TOPAS\DS_carbon_MultiEnergyLocal\LUTResult\';
doubleS.d = [1:160]*3 - 3/2;%[1:240]*2 - 1;

doubleS.Ions = {'protons', 'He', 'Li', 'Be', 'B', 'C'};

for k=1:size(doubleS.Ions,2)
   doubleS.EParam(k).EMax = 250.1;
   doubleS.EParam(k).EMin = 0.1;
   doubleS.EParam(k).nEBins = 250;
   
   doubleS.EParam(k).EdMax = 12;
   doubleS.EParam(k).EdMin = 0;
   doubleS.EParam(k).nEdBins = 240;
end

doubleS.computeBinning();

doubleS.rn = 0.25;
doubleS.RN = 3.9;
doubleS.Beta_Tissue = 0.05;

doubleS.analyzeRawData(0);

doubleS.LUT = 'D:\matRad_gitHubRemo_MKM_TOPAS\MKM\LUTs_rn25_RN_39_Beta005_500MeV';
doubleS.computeZD();
FirstLUT = doubleS.zD;

doubleS.LUT = 'D:\matRad_gitHubRemo_MKM_TOPAS\MKM\LUTs_rn35_RN31_Beta005_250MeV';
doubleS.computeZD();
SecondLUT = doubleS.zD;

doubleS.LUT = 'D:\matRad_gitHubRemo_MKM_TOPAS\MKM\LUTs_rn27_RN39_Beta005_300MeV';

doubleS.computeZD();
ThirdLUT = doubleS.zD;
%% Plot
figure;
plot(doubleS.d, FirstLUT, '.-');
hold on;

plot(doubleS.d, SecondLUT, '.-');
plot(doubleS.d, ThirdLUT, '.-');

legend('rn = 0.25', 'rn = 0.35');
grid on;
%% Base Data gen
doubleS = DS_MultiIon();

doubleS.wDir = 'D:\matRad_gitHubRemo_MKM_TOPAS\DS_Carbons_MultiEnergyCluster_Long\LUTResultLong';
doubleS.d = [1:240]*2 - 1;%[1:240]*2 - 1;
doubleS.LUT = 'D:\matRad_gitHubRemo_MKM_TOPAS\MKM\LUTs_rn25_RN_39_Beta005_500MeV';

doubleS.Ions = {'protons', 'He', 'Li', 'Be', 'B', 'C'};


for k=1:size(doubleS.Ions,2)
   doubleS.EParam(k).EMax = 300.0;
   doubleS.EParam(k).EMin = 0.1;
   doubleS.EParam(k).nEBins = 300;
   
   doubleS.EParam(k).EdMax = 12;
   doubleS.EParam(k).EdMin = 0;
   doubleS.EParam(k).nEdBins = 240;
end

doubleS.computeBinning();

for k=1:28
   doubleS.analyzeRawData(k-1);
   doubleS.computeZD();

   %%%%% !!!!!! %%%% 3 mm offset inserted
   Spe(k).d = doubleS.d;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Spe(k).zMix  = doubleS.zD;
   Spe(k).DS.depths = doubleS.d;
   Spe(k).DS.EBinning = doubleS.EParam(1).E;
   Spe(k).DS.EdBinning = doubleS.EParam(1).Ed;
   Spe(k).DS.Phi       = doubleS.CumulativeEdPhi;
   display(['Computing zD for Energy: ', num2str(k+1)]);


end
%

figure;
for k=1:10
   plot(Spe(k).d, Spe(k).zMix, '.-');

   hold on;

end


%% Plot LUT
figure;
for k=1:6
   lut = load(strcat(doubleS.LUT, filesep, 'LUT_Z_', num2str(k), '.mat'));
   ELUT = lut.E;
   zDLUT = lut.zD;
   semilogx(ELUT,zDLUT, '.-');
   hold on;
end

legend(doubleS.Ions);
%% Insert machine data
[~,Ex] = intersect([machine.data(:).energy],stf.ray(1).energy);
for k=1:size(Ex,1)
     
       machine.data(Ex(k)).Spe = Spe(k);
end

%% MKM vs LEM

M = load(['carbon_GenericMKM_DS_rn25.mat']);
M = M.machine;


EIdx = 31;

depths = M.data(EIdx).depths;
LEM_Param = matRad_bioModel('carbon', 'RBExD', 'LEM');
LEM_Profile = LEM_Param.calcLQParameter(depths,M.data(EIdx),ones(size(depths)),0.1*ones(size(depths)),0.05*ones(size(depths)),2*ones(size(depths)));

stf = matRad_computeSSD(stf,ct);
BAMStoIsoDist = 2000;
nozzleToSkin = ((stf.ray(1).SSD + BAMStoIsoDist) - machine.meta.SAD);
dR = 0.0011 * (nozzleToSkin);

depths_1 = M.data(EIdx).depths;% + dR;
MKM_Param = matRad_bioModel('carbon', 'RBExD', 'MKM');
MKM_Profile = MKM_Param.calcLQParameter(depths_1,M.data(EIdx),ones(size(depths_1)),0.1*ones(size(depths_1)),0.05*ones(size(depths)),2*ones(size(depths_1)));
figure;
plot(depths, LEM_Profile./LEM_Profile(2), '.-');
hold on;
plot(depths_1, MKM_Profile./MKM_Profile(2), '.-');
grid on;
legend('LEM', 'MKM');