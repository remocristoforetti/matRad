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

%% Dose calculation
% 
% load('D:\matRad_gitHubRemo_MKM_TOPAS\MKM\Spectrum_CarbonsMultiEnergy.mat');
% [~,Ex] = intersect([machine.data(:).energy],stf.ray(1).energy);
% for k=1:size(Ex,1)
%      
%        machine.data(Ex(k)).Spe = Spe(k);
% end
   
if (strcmp(modelName, 'MKM') && ~strcmp(pln.machine,'GenericMKM'))
   try
      load([pln.radiationMode, '_', pln.machine, 'MKM']);
      pln.machine = 'GenericMKM';
   catch
      matRad_cfg.dispWarning('Unable to find machine data');
   end
end

AvailableDataIdx = find(arrayfun(@(x) ~isempty(x.Spe), [machine.data]));
if ~size(intersect([machine.data(AvailableDataIdx).energy],unique([stf.ray(:).energy])),2) == size(unique(stf.ray(1).energy),2)
   matRad_cfg.dispWarning('Energies are not correct');
end

pln.machine = 'GenericMKM';
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
plot(x, Profile_LEM, '.-');
hold on;
plot(x, Profile_MKM, '.-');
plot(x, Profile_pd, '.-');
ylabel('Dose [Gy]');
yyaxis right;
plot(x, Profile_MKM_LET, '.-');
grid on;
ylabel('KeV/\mu');
xlabel('Depth [mm]');
legend('LEM', 'MKM', 'physical Dose', 'LET');
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

topas_cfg.writeAllFiles(ctT,cst,pln,stf,machine,w);
WriteLUT_Multition('X:\matRad_dev_varRBE_robOpt_TopasMKM\MKM\LUTs', topas_cfg);

AddSubscorersSingleIons = 0;
if AddSubscorersSingleIons
    Ions = [1 2 3 4 5 6];
    WriteSubScorers(topas_cfg,Ions);
end
%% ReadOut TopasData
filenames = dir('D:\matRad_gitHubRemo_MKM_TOPAS\MKM\MKM_SimulationData_Carbon_Tot\*.dcm');
zMix = [];
d = [1:480] - 0.5;%[1:480] - 0.5; %[1:160]*3 -3/2;%-0.5;
for k=1:size(filenames)
    infos = dicominfo(strcat(filenames(1).folder, filesep, filenames(k).name));
    zMix = [zMix, double(dicomread(strcat(filenames(1).folder, filesep, filenames(k).name)))*infos.DoseGridScaling];
end

% filenames = dir('C:\Users\Remo\Desktop\SimulationOutput\Carbon_LUT_Tot_Weight\*.dcm');
% zMix_1 = [];
% d = [1:480] - 0.5;%[1:480] - 0.5; %[1:160]*3 -3/2;%-0.5;
% for k=1:size(filenames)
%     infos = dicominfo(strcat(filenames(1).folder, filesep, filenames(k).name));
%     zMix_1 = [zMix_1, double(dicomread(strcat(filenames(1).folder, filesep, filenames(k).name)))*infos.DoseGridScaling];
% end


toplot = [1:5:28];

figure;
plot(d, zMix(:,toplot), '.-');
% hold on;
% plot(d, zMix_1(:,toplot), '.-');
grid on;
xlim([0 200]);

% for k = 1:size(zMix,2)
%     Spe(k).d = d;
%     Spe(k).zMix = zMix(:,k);
% end
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

%% Plot Fragments

DataFiles = dir('D:\matRad_gitHubRemo_MKM_TOPAS\Fragment results\*.dcm');
lege = [];
d = [1:480] - 0.5;
figure;
for k=1:size(DataFiles,1)
    Filename = strcat(DataFiles(1).folder, filesep, DataFiles(k).name);
    Data = dicomread(Filename);
    DataInfo = dicominfo(Filename);
    
    if ~strcmp(DataFiles(k).name, 'ZMixOutput.dcm')
        IonName = erase(DataFiles(k).name, 'ZMixOutput_');
        IonName = erase(IonName, '.dcm');
        lege = [lege, {IonName}];
    else
        lege = [lege, {'Total ZMix'}];
    end

    plot(d, double(Data)*DataInfo.DoseGridScaling, '.-');
    hold on;

end
xlim([0 200]);
grid on;
legend(lege);

%% Read Double Spectrum
v = dicomread('D:\matRad_gitHubRemo_MKM_TOPAS\Ds.dcm');
infos = dicominfo('D:\matRad_gitHubRemo_MKM_TOPAS\Ds.dcm');
v = double(squeeze(v));%*infos.DoseGridScaling;
m= permute(v, [2,3,1]);

toplot = [22:1:32];
for k=toplot
   figure;
   imagesc(double(m(2:end,2:end,k))./sum(m, [1,2,3]));
   title('Z = ',num2str(k*3 +3/2));
   currXtick = get(gca, 'XTick');
   currYtick = get(gca, 'YTick');

   set(gca,'XTick',currXtick, 'XTickLabel', linspace(currXtick(1)*150/(currXtick(end)+10),150,size(currXtick,2)));

   set(gca,'YTick', currYtick, 'YTickLabel', linspace(currYtick(1)*5/(currYtick(end)),5.0,size(currYtick,2)));
   xlabel('E [MeV]');
   ylabel('Edep [MeV]');

   colorbar;
   
   
end
