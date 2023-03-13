matRad_rc;
matRad_cfg = MatRad_Config.instance();
%% Setup dummy Phantom
resolution.x = 1; %mm
resolution.y = 1; %mm
resolution.z = 1; %mm

cubeRealDim = [300 100 100]; %mm
cubeDim  = ceil(cubeRealDim./[resolution.y resolution.x resolution.z]); %In Voxels

ct.cube = {ones(cubeDim)};
ct.resolution = resolution;
ct.cubeDim = cubeDim;
ct.numOfCtScen = 1;

ct.cubeHU = {zeros(cubeDim)}; %All water

ct.hlut = [1 0; ...
           0 -1024 ];
ct.ctGrid.resolution = ct.resolution;
%%%%%%% Example %%%%%%%%%%%%%

% % % resolution.x = 1; %mm
% % % resolution.y = 0.2; %mm
% % % resolution.z = 1; %mm
% % % 
% % % cubeRealDim = [150 100 100]; %mm
% % % cubeDim  = ceil(cubeRealDim./[resolution.y resolution.x resolution.z]); %In Voxels
% % % 
% % % ct.cube = {ones(cubeDim)};
% % % ct.resolution = resolution;
% % % ct.cubeDim = cubeDim;
% % % ct.numOfCtScen = 1;
% % % 
% % % ct.cubeHU = {zeros(cubeDim)}; %All water
% % % 
% % % ct.hlut = [1 0; ...
% % %            0 -1024 ];
% % % ct.ctGrid.resolution = ct.resolution;
%% Setup

%pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = 0;
pln.propStf.couchAngles   = 0;
pln.propStf.bixelWidth    = 3;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
%pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

pln.radiationMode   = 'protons';
pln.propDoseCalc.calcLET = 1;

%% Dummy machine
% machine.meta.radiationMode      = pln.radiationMode;
% machine.meta.dataType           = 'singleGauss';
% machine.meta.created_on         = 'never';
% machine.meta.created_by         = 'nobody';
% machine.meta.description        = 'dummy machine';
% machine.meta.SAD                = 10000;
% machine.meta.BAMStoIsoDist       = 2000;
% machine.meta.LUT_bxWidthminFWHM = [1 Inf,...
%                                    3 3];
% machine.meta.machine            = 'dummy';

%Retrive some infos from Generic machine
%Generic_machine = load('carbon_Generic.mat');
%HITmachine = load('carbon_HITgantry');
topasMachine = load('protons_generic_TOPAS_mcSquare.mat');
%genericMachine = load('protons_Generic.mat');
%machine = Generic_machine.machine;
%machine = HITmachine.machine;
machine = topasMachine.machine;
%machine = genericMachine.machine;
%% Dummy STF

%For the time being, keep the same energies saved in the Generic machine
%E = [machine.data(31:37).energy];
E = [machine.data(19:41).energy];


stf.gantryAngle   = pln.propStf.gantryAngles;
stf.couchAngle    = pln.propStf.couchAngles;
stf.bixelWidth    = pln.propStf.bixelWidth;
stf.radiationMode = pln.radiationMode;
stf.SAD           = machine.meta.SAD;
% stf.isoCenter = [0.5*ct.resolution.x*(ct.cubeDim(2)), ...
%                   ct.resolution.y*(ct.cubeDim(1) +1), ...
%                   0.5*ct.resolution.z*(ct.cubeDim(3))];
stf.isoCenter = [0.5*ct.resolution.x*(ct.cubeDim(2)), ...
                  0, ...
                  0.5*ct.resolution.z*(ct.cubeDim(3))];

stf.numOfRays = 1;

   ray.rayPos_bev       = [0 0 0]; 
   ray.targetPoint_bev  = [0 stf.SAD 0];
   ray.rayPos           = [0 0 0];
   ray.targetPoint      = [0 stf.SAD 0];
   ray.energy           = E;
   ray.rangeShifter.ID  = 0;
   ray.rangeShifter.eqThickness = 0;
   ray.rangeShifter.sourceRashiDistance = 0;
   ray.focusIx = ones(1,numel(E));
   ray.numParticlesPerMU = 1000000*ones(1,34);
   ray.minMU = zeros(size(ray.numParticlesPerMU));
   ray.MaxMU = Inf(size(ray.numParticlesPerMU));

   stf.ray = ray;
stf.sourcePoint_bev        = [0 -stf.SAD 0];
stf.sourcePoint            = [0 -stf.SAD 0];

stf.numOfBixelsPerRay      = numel(stf.ray.energy);
stf.longitudinalSpotSpacig = 2;

stf.totalNumOfBixels        = sum(stf.numOfBixelsPerRay);

switch pln.radiationMode
   case 'protons'
      Ions = {'protons'};
   case 'carbon'
      Ions = {'protons', 'He', 'Li', 'Be', 'B', 'C'};
end
for k=1:size(Ions,2)
   EParam(k).EMax = 250;
   EParam(k).EMin = 0;
   EParam(k).nEBins = 500;

   EParam(k).EdMax = 12;
   EParam(k).EdMin = 0;
   EParam(k).nEdBins = 240;
end
% for k=1:size(Ions,2)
%    EParam(k).EMax = 230.1;
%    EParam(k).EMin = 0.1;
%    EParam(k).nEBins = 460;
% 
%    EParam(k).EdMax = 12;
%    EParam(k).EdMin = 0;
%    EParam(k).nEdBins = 120;
% end
if any([EParam.nEdBins] == 1)
   matRad_cfg.dispWarning('Ed binning deactivated');
end
%% TOPAS SIMULATION

topas_cfg = matRad_TopasConfig;
topas_cfg.worldMaterial = 'Vacuum';
folderName = 'BaseData_proton_spectra/SOBP';
if ~exist(strcat(topas_cfg.workingDir, folderName), 'dir')
   mkdir(strcat(topas_cfg.workingDir, folderName));

end
topas_cfg.workingDir = strcat(topas_cfg.workingDir, folderName);

topas_cfg.label = 'BaseData';
topas_cfg.numOfRuns = 1;
topas_cfg.numThreads = 0;
numHistoriesPerBixel = 100000;

topas_cfg.numHistories = ceil(numHistoriesPerBixel*size(stf.ray.energy,2));

includePDD = true;
includeDS  = false;
includeSP  = true;

includeEdEventSpectrum = true;
includeProtonLET = true;
nPhantoms = 3;
if exist('nPhantoms', 'var')
    highResWindowWidtgh = 30;  %in mm

    distalPhantomsWidth = 30;  %in mm
    resX = 5*ones(1,3);
    resY = [0.1 0.1 0.1];
    resZ = 5*ones(1,3);
    proximalPhantomWidth = ceil(machine.data(31).peakPos - highResWindowWidtgh/2);
    dimensionsDepth = [proximalPhantomWidth,highResWindowWidtgh,distalPhantomsWidth];
    positionDepth   = cumsum(dimensionsDepth);
    cubeDimensions  = [sum(dimensionsDepth) 100 100]; %In mm
    for k=1:nPhantoms
      phantoms(k).name         = ['Phantom' num2str(k)];
      phantoms(k).dimension    = [dimensionsDepth(k), cubeDimensions(2), cubeDimensions(3)];
      phantoms(k).resolution   = [resY(k), resX(k), resZ(k)];
    
      phantoms(k).positionDepth = [0, positionDepth];
      phantoms(k).EParam        = EParam;

    end
    writeScoringTreeDirectory(phantoms,Ions,topas_cfg,EParam,includePDD,includeDS,includeSP,includeEdEventSpectrum,includeProtonLET);
end


%Define highRes phantoms
% % % % % % nPhantoms = 3;
% % % % % % 
% % % % % % if exist('nPhantoms', 'var')
% % % % % % 
% % % % % %    cubeDimensions  = [150 100 100]; %In mm
% % % % % %    dimensionsDepth = [75 55 20];    %In mm, this is the phantom dimension in depth
% % % % % %    positionDepth = cumsum(dimensionsDepth);
% % % % % %    resX = 5*ones(1,3);
% % % % % %    resY = [1 0.1 1];
% % % % % %    resZ = 5*ones(1,3);
% % % % % % 
% % % % % %    for k=1:nPhantoms
% % % % % %       phantoms(k).name         = ['Phantom' num2str(k)];
% % % % % %       phantoms(k).dimension    = [dimensionsDepth(k), cubeDimensions(2), cubeDimensions(3)];
% % % % % %       phantoms(k).resolution   = [resY(k), resX(k), resZ(k)];
% % % % % % 
% % % % % %       phantoms(k).positionDepth = [0, positionDepth];
% % % % % %       phantoms(k).EParam        = EParam;
% % % % % %    end
% % % % % % 
% % % % % %    writeScoringTreeDirectory(phantoms, Ions,topas_cfg, EParam);
% % % % % % end
   
%writeDSscorers(EParam, Ions,topas_cfg);

%topas_cfg.scorer.LET             = pln.propDoseCalc.calcLET;


topas_cfg.scorer.outputType      = 'DICOM';
topas_cfg.scorer.reportQuantity  = {'Sum'};
topas_cfg.scorer.volume          = true;
topas_cfg.scorer.doseToMedium    = false; %Turn down and add parallel scorers for higher resolution
topas_cfg.scorer.filename        = 'constructor';

w = (1/size(stf.ray.energy,2))*ones(size(stf.ray.energy,2),1);

topas_cfg.writeAllFiles(ct,0,pln,stf,machine,w);

% if exist('nPhantoms', 'var')
% 
%    for k=1:nPhantoms
%        if ~exist(strcat(topas_cfg.workingDir,filesep, 'Output\PDD',filesep, phantoms(k).name), 'dir')
%          mkdir(strcat(topas_cfg.workingDir,filesep, 'Output\PDD',filesep, phantoms(k).name));
%       end
%    end
% 
%    for k=1:6
% 
%       for m=1:nPhantoms
%           if ~exist(strcat(topas_cfg.workingDir,filesep, 'Output\DS',filesep, phantoms(m).name,filesep, strcat('Ion_',num2str(k))), 'dir')
%             mkdir(strcat(topas_cfg.workingDir,filesep, 'Output\DS',filesep, phantoms(m).name,filesep, strcat('Ion_',num2str(k))));
%          end
%       end
%       
%    end
% end

%% HighRes phantom data acquisistion
analisis = highResPhantomAnalisis();

analisis.nPhantoms = nPhantoms;
analisis.wDir     = 'C:\Users\r408i\Desktop\r408i_data\BaseData_protons_spectra\SOBP'; %topas_cfg.workingDir;
analisis.readPDD   = true;
analisis.readSpectra = true;
analisis.readDSSpectra = false;
analisis.readSPeDSpectra = true;
analisis.readProtonLET = true;
analisis.ions = {'protons'};
analisis.integrationSurface = [1:20];
%Do this after ions have been defined
analisis.phantoms  = phantoms;



for k=1:size(E,2)
    analisis.run       = k-1;
    analisis.importRawData();
    DS(k).Phi       = analisis.Phi;
    DS(k).PDD       = analisis.PDD;
    DS(k).depths    = analisis.depths;
    DS(k).edPhi     = analisis.edPhi;

    DS(k).eBinning  = analisis.phantoms(1).EParam(1).E;
    DS(k).ProtonLET = analisis.ProtonLET;
end
%% Compare baseData and new PDD
figure;
for k=1:10:size(E,2)
    plot(machine.data(19 + k -1).depths,machine.data(19 + k -1).Z./max(machine.data(19  + k - 1).Z()), '--');
    hold on;
    plot(analisis.depths, DS(k).PDD./max(DS(k).PDD()), '.-');
    grid on;
end
%% Load into baseData
newMachine = machine;
[~,eIdx] = intersect([machine.data.energy], E);
for k=1:size(eIdx,1)
    newMachine.data(eIdx(k)).DS.Phi = DS(k).Phi;
    newMachine.data(eIdx(k)).DS.edPhi = DS(k).edPhi;
    newMachine.data(eIdx(k)).DS.eBinning = DS(k).eBinning; %This is valid only when resolutions are all the same. Need to handel different resolutions.
    newMachine.data(eIdx(k)).DS.depths = DS(k).depths;
    newMachine.data(eIdx(k)).DS.PDD = DS(k).PDD;
    newMachine.data(eIdx(k)).LET = interp1(DS(k).depths, DS(k).ProtonLET, machine.data(eIdx(k)).depths);
end
% newMachine = machine;
% [~,eIdx] = intersect([machine.data.energy], E);
% for k=1:size(eIdx,1)
%     newMachine.data(eIdx(k)).DS.Phi = analisis.Phi;
%     newMachine.data(eIdx(k)).DS.edPhi = analisis.edPhi;
%     newMachine.data(eIdx(k)).DS.eBinning = analisis.phantoms(1).EParam(1).E; %This is valid only when resolutions are all the same. Need to handel different resolutions.
%     newMachine.data(eIdx(k)).DS.depths = analisis.depths;
%     newMachine.data(eIdx(k)).DS.PDD = analisis.PDD;
%     newMachine.data(eIdx(k)).LET = interp1(analisis.depths, analisis.ProtonLET, machine.data(eIdx(k)).depths);
% end
% oldMachine = machine;
% machine = newMachine;
% save(['protons_newGenericTOPAS.mat'], 'machine');
% machine = oldMachine;
%% Spectra analiisis
DSPhi = full(analisis.DSPhi{6});
Phi   = full(analisis.Phi{6});

depthProfile = 910;
plEnergies = analisis.phantoms(1).EParam(1).E;

figure;
plot(plEnergies, DSPhi(:,depthProfile), '.-');
hold on;
plot(plEnergies, Phi(:,depthProfile), '.-');
grid on;
yyaxis 'right'
plot(plEnergies, DSPhi(:,depthProfile)./Phi(:,depthProfile), '.-');
ylim([1 50]);
legend('event fluence', 'particle fluence');

% ratios = [];
% for k=1:size(analisis.depths,2)
%     maxDSPhi = max(DSPhi(:,k));
%     maxPhi = max(Phi(:,k));
% 
%     ratios = [ratios, maxDSPhi/maxPhi];
% end
% 
% figure;
% plot(analisis.depths, ratios, '.-');
% xlim([0,100]);
%% E-R analisis
baseDataEnergies = [machine.data(:).energy];

baseDataRanges = [];
%%%% MCemittance r80 sampling
for k=1:size(baseDataEnergies,2)
   newDepths = linspace(0,machine.data(k).depths(end),numel(machine.data(k).depths) * 100);
   newDose   = interp1(machine.data(k).depths, machine.data(k).Z, newDepths, 'spline');

   [maxV, maxI] = max(newDose);
   [~, r80ind] = min(abs(newDose(maxI:end) - 0.8 * maxV));
   r80ind = r80ind - 1;


   
   baseDataRanges(k) = interp1(newDose(maxI + r80ind - 1:maxI + r80ind + 1), ...
                  newDepths(maxI + r80ind - 1:maxI + r80ind + 1), 0.8 * maxV);

end

meanEnergy = @(x) 11.39 * x^0.628 + 11.24;
fitEnergies = arrayfun(@(r80) meanEnergy(r80), baseDataRanges);

fitRanges = interp1(fitEnergies, baseDataRanges, baseDataEnergies, 'spline');

figure;
subplot(2,1,1);
plot(baseDataEnergies, baseDataRanges, '.-');
hold on;
plot(baseDataEnergies, fitRanges, '.-');
grid on;
xlabel('Energy [MeV]');

ylabel('Range [mm]');
legend('baseData', 'fit');
subplot(2,1,2);

plot(baseDataEnergies, baseDataRanges - fitRanges, '.-');
% hold on;
% plot(baseDataEnergies, baseDataRanges - [machine.data(:).range], '.-');

grid on;
xlabel('Energy [MeV]');
ylabel('baseData - fit [mm]');
%% E-R Repeat fit
for k=1:size(baseDataEnergies,2)
   newDepths = linspace(0,machine.data(k).depths(end),numel(machine.data(k).depths) * 100);
   newDose   = interp1(machine.data(k).depths, machine.data(k).Z, newDepths, 'spline');

   [maxV, maxI] = max(newDose);
   [~, r80ind] = min(abs(newDose(maxI:end) - 0.8 * maxV));
   r80ind = r80ind - 1;


   
   baseDataRanges(k) = interp1(newDose(maxI + r80ind - 1:maxI + r80ind + 1), ...
                  newDepths(maxI + r80ind - 1:maxI + r80ind + 1), 0.8 * maxV);

end
%baseDataRanges = baseDataRanges + machine.data(1).offset;
F = fittype('a + b*x^c + d*x^e', 'coeff', {'a', 'b', 'c', 'd', 'e'});
outFit = fit(baseDataRanges', baseDataEnergies', F, 'StartPoint', [11.24, 11.39, 0.628, 1,1]);

newFitEnergies = outFit(baseDataRanges);
newFitRanges = interp1(newFitEnergies, baseDataRanges, baseDataEnergies, 'spline');

figure;
subplot(2,1,1);

plot(baseDataEnergies, baseDataRanges, '.-');
hold on;
plot(baseDataEnergies, newFitRanges, '.-');
grid on;
xlabel('Energy [MeV]');

ylabel('Range [mm]');
legend('baseData', 'new fit');
subplot(2,1,2);
plot(baseDataEnergies, baseDataRanges - newFitRanges, '.-');
grid on;
xlabel('Energy [MeV]');
ylabel('baseData - fit [mm]');

%% E-R: TOPAS comparison
eIdx = [1:10:size(machine.data,2)];
eIdx = eIdx(1:end-2);
energies = [machine.data(eIdx).energy];
addOffset = false;
wDir = 'C:\Users\r408i\Desktop\r408i_data\carbonBaseDataChecks\generic\path2';
name = 'score_BaseData_field1_run1_physicalDose_WATER_Run_';

x = [1:1500]*0.2 - 0.1;
integration = [40:60];
%Load vaccum profiles
PDD = [];
topasR80_HIT = [];
for k=0:size(eIdx,2)-1
   fileName = strcat(wDir, filesep,name, string(compose('%04d', k)));
   data = squeeze(double(dicomread(fileName)));
   data = permute(data, [2,3,1]);

   profileData = squeeze(sum(data(integration,integration,:),[1,2]));
   
   [maxV, maxI] = max(profileData);
   [~, r80ind] = min(abs(profileData(maxI:end) - 0.8 * maxV));
   r80ind = r80ind - 1;
   topasR80_HIT(k+1) = interp1(profileData(maxI + r80ind - 1:maxI + r80ind + 1), ...
             x(maxI + r80ind - 1:maxI + r80ind + 1), 0.8 * maxV);
    if addOffset
      topasR80_HIT(k+1) = topasR80_HIT(k+1) + machine.data(k+1).offset;
   end
   PDD(:,k+1) = profileData;
end

figure;
lege = [];
for k=1:size(eIdx,2)
   plot(x, PDD(:,k)./PDD(1,k), '.-');
   hold on;
   plot(machine.data(eIdx(k)).depths, machine.data(eIdx(k)).Z./machine.data(eIdx(k)).Z(1), '--');
   lege = [lege, {[num2str(energies(k)), 'MeV Topas']}, {'Nominal'}];
end
grid on;
xlabel('Depth [mm]');
ylabel('Dose');
legend(lege);
title('HIT - Nominal energies');
%% E-R TOAPS, Comparison
figure;
subplot(2,1,1);
plot([machine.data(eIdx).energy],baseDataRanges(eIdx), 'o-');
hold on;

plot([machine.data(eIdx).energy],topasR80_HIT, 'o-');

xlabel('Energy [MeV]');
ylabel('Range [mm]');
grid on;
legend('baseData', 'TOPAS Vacuum');
title('E-R relation, HIT baseData');
subplot(2,1,2);
plot([machine.data(eIdx).energy],baseDataRanges(eIdx)-topasR80_HIT, 'o-');
% hold on;
% plot([machine.data(eIdx).energy],baseDataRanges(eIdx)-machine.data(1).offset-topasR80_HIT, 'o-');

xlabel('Energy [MeV]');
ylabel('baseData - topas [mm]');
grid on;

%% Repeat fit for generic base data

%% Base data fluence analisis
resultDir = 'D:\matRad_gitHubRemo_MKM_TOPAS\topas\Results\BaseData_Carbon7\';

doubleS = DS_MultiIon();
doubleS.wDir = resultDir;
doubleS.Ions = {'protons', 'He', 'Li', 'Be', 'B','C'};
doubleS.d = [1:750]*ct.resolution.y - ct.resolution.y/2;
doubleS.EParam = EParam;
doubleS.computeBinning();
for k=1:size(E,2)
   display(['Computing zD for Energy: ', num2str(k+1)]);
   doubleS.analyzeRawData(k-1);

   doubleS.readDoseProfile(k-1,1);
   PDD(:,k) = doubleS.PDD;
end

figure;
plot(doubleS.d, PDD(:,7), '.-');
%% Add DS to BaseData
%For the time being, just load the DS into BaseData and move to MKM
%calculation
resultDir = 'D:\matRad_gitHubRemo_MKM_TOPAS\topas\MCrun\BaseData';

doubleS = DS_MultiIon();

doubleS.wDir = resultDir;
doubleS.d = linspace(0,ct.resolution.y*(ct.cubeDim(1)),cubeDim(1)) + 0.5*ct.resolution.y;

doubleS.LUT = 'D:\matRad_gitHubRemo_MKM_TOPAS\MKM\LUTs_rn25_RN_39_Beta005_500MeV';

doubleS.Ions = {'protons', 'He', 'C'};

for k=1:size(doubleS.Ions,2)
   doubleS.EParam(k).EMax = 250.1;
   doubleS.EParam(k).EMin = 0.1;
   doubleS.EParam(k).nEBins = 250;
   
   doubleS.EParam(k).EdMax = 12;
   doubleS.EParam(k).EdMin = 0;
   doubleS.EParam(k).nEdBins = 240;
end

doubleS.computeBinning();

for k=1:size(E,2)
   doubleS.analyzeRawData(k-1);
   %doubleS.computeZD();
   %DS.zD = doubleS.zD;
   display(['Computing zD for Energy: ', num2str(k+1)]);


end
[~,MIdx] = intersect([machine.data(:).energy],[BaseData.energy]);
for k=1:size(E,2)
   machine.data(MIdx).DS.IonCharges = doubleS.IonNumber;
   machine.data(MIdx).DS.EBins      = doubleS.EParam.E;
   machine.data(MIdx).DS.EdBins     = doubleS.EParam.Ed;
   machine.data(MIdx).DS.Fluence    = doubleS.SparsePhi;
end

%
% figure;
% for k=1:10
%    plot(Spe(k).d, Spe(k).zMix, '.-');
% 
%    hold on;
% end



%% Data Import and Analysis
%Complete BaseData creation should go here
ResultDir = 'D:\matRad_gitHubRemo_MKM_TOPAS\topas\MCrun\BaseData';

DoseScorer = 'score_BaseData_field1_run1_physicalDose_ES.dcm';
BaseData = [];
for Eidx = 1:size(E,2)
   DoseData = dicomread([ResultDir,filesep,DoseScorer]);
   DoseInfo = dicominfo([ResultDir,filesep,DoseScorer]);
   DoseData = squeeze(DoseData);
   DoseData = double(DoseData)*DoseInfo.DoseGridScaling;
   %DoseData = permute(DoseData, [2,3,1]);
   resolution.x = ct.resolution.y;
   resolution.y = ct.resolution.x;
   resolution.z = ct.resolution.z;
   fData = matRad_fitBaseData(DoseData(:,:,:),resolution,E(Eidx));
   BaseData(Eidx).energy   = fData.energy;
   BaseData(Eidx).range    = fData.range;
   
   BaseData(Eidx).depths   = fData.depths;
   BaseData(Eidx).Z        = fData.Z;
   
   BaseData(Eidx).peakPos  = fData.peakPos;
   BaseData(Eidx).offset   = fData.offset;
   BaseData(Eidx).sigma1   = fData.sigma1;
   BaseData(Eidx).sigma2   = fData.sigma2;

%    BaseData(Eidx).energy = E(Eidx);

   
   BaseData(Eidx).MCdepths = linspace(0,ct.resolution.y*(ct.cubeDim(1)),cubeDim(1))' + 0.5*ct.resolution.y;
   cf = 1 / 1.6021766208e-02 * ct.resolution.y * ct.resolution.z;
   BaseData(Eidx).MCZ = squeeze(sum(DoseData, [2,3]))*cf;
%    
%    
%    BaseData(Eidx).offset = 0;
   
end

%Compute and insert DS

doubleS = DS_MultiIon();
doubleS.wDir = resultDir;
doubleS.d = BaseData(1).depths;
doubleS.LUT = 'D:\matRad_gitHubRemo_MKM_TOPAS\MKM\LUTs_rn25_RN_39_Beta005_500MeV';
doubleS.Ions = {'protons', 'He', 'Li', 'Be', 'B', 'C'};

for k=1:size(doubleS.Ions,2)
   doubleS.EParam(k).EMax = 250.1;
   doubleS.EParam(k).EMin = 0.1;
   doubleS.EParam(k).nEBins = 250.1;

   doubleS.EParam(k).EdMax = 12;
   doubleS.EParam(k).EdMin = 0;
   doubleS.EParam(k).nEdBins = 240;

end
doubleS.computeBinning();

[~,MIdx] = intersect([machine.data(:).energy],[BaseData.energy]);

for Eidx =1:size(E,2)
   display(['Computing Fluence for Energy: ', num2str(Eidx)]);
   doubleS.analyzeRawData(Eidx-1);
   machine.data(MIdx(Eidx)).DS.IonCharges = doubleS.IonNumber;
   machine.data(MIdx(Eidx)).DS.EBins      = doubleS.EParam.E;
   machine.data(MIdx(Eidx)).DS.EdBins     = doubleS.EParam.Ed;
   machine.data(MIdx(Eidx)).DS.Fluence    = doubleS.SparsePhi;

end

%% Plot
%[~,MIdx] = intersect([machine.data(:).energy],[BaseData.energy]);
MIdx = 31;
figure;

plot(machine.data(MIdx).depths, machine.data(MIdx).Z./machine.data(MIdx).Z(1), '.-');
hold on;
plot(BaseData.MCdepths-0.5, BaseData.MCZ./BaseData.MCZ(1), '.-');
grid on;
xlim([70 100]);
%plot(BaseData.depths - BaseData.offset, BaseData.Z./BaseData.Z(1), '.-');
%plot(BaseData.depths - BaseData.offset, BaseData.Z, '.-');

%plot(machine.data(MIdx).depths,machine.data(MIdx).Z, '.-');
%
%plot(machine.data(MIdx).depths,machine.data(MIdx).Z./machine.data(MIdx).Z(1), '.-');
%hold on;
%plot(fData.depths- BaseData.offset, fData.Z./fData.Z(1), '.-');

%plot(fData.depths- BaseData.offset, fData.Z, '.-');
%grid on;
%xlim([-1 110]);
%% Plot PDD vs BaseData
%PDD = 'D:\matRad_gitHubRemo_MKM_TOPAS\PDD_protons_vacuum.dcm';
%PDD = double(dicomread(PDD));
%PDD_phantom = 'D:\matRad_gitHubRemo_MKM_TOPAS\PDD_protons_vacuum_phantom.dcm';
integr = [80:80];
PDD_vacuum_nominal27 = 'D:\matRad_gitHubRemo_MKM_TOPAS\PDD_protons_vacuum_phantom.dcm';
PDD_vacuum_nominal27 = squeeze(double(dicomread(PDD_vacuum_nominal27)));
PDD_vacuum_nominal27 = permute(PDD_vacuum_nominal27, [2,3,1]);
PDD_vacuum_nominal27 = squeeze(sum(PDD_vacuum_nominal27(integr,integr,:), [1,2]));

PDD_air_nominal27 = 'D:\matRad_gitHubRemo_MKM_TOPAS\PDD_protons_air_nominal27.dcm';
PDD_air_nominal27 = squeeze(double(dicomread(PDD_air_nominal27)));
PDD_air_nominal27 = permute(PDD_air_nominal27, [2,3,1]);
PDD_air_nominal27 = squeeze(sum(PDD_air_nominal27(integr,integr,:), [1,2]));


PDD_vacuum_MC27 = 'D:\matRad_gitHubRemo_MKM_TOPAS\PDD_protons_vacuum_MCemittance27.dcm';
PDD_vacuum_MC27 = squeeze(double(dicomread(PDD_vacuum_MC27)));
PDD_vacuum_MC27 = permute(PDD_vacuum_MC27, [2,3,1]);
PDD_vacuum_MC27 = squeeze(sum(PDD_vacuum_MC27(integr,integr,:), [1,2]));

PDD_air_MC27 = 'D:\matRad_gitHubRemo_MKM_TOPAS\PDD_protons_air_MCemitt27.dcm';
PDD_air_MC27 = squeeze(double(dicomread(PDD_air_MC27)));
PDD_air_MC27 = permute(PDD_air_MC27, [2,3,1]);
PDD_air_MC27 = squeeze(sum(PDD_air_MC27(integr,integr,:), [1,2]));

%% plot

load('protons_Generic.mat');
PDD_basedata = machine.data(27).Z;

leg = [];
figure;
% plot([1:480] -0.5, PDD_vacuum_nominal27./max(PDD_vacuum_nominal27), '.-');
% leg = [leg, {'Vacuum nominal'}];
hold on;

plot([1:480] -0.5, PDD_air_nominal27./max(PDD_air_nominal27), '.-');
leg = [leg, {'Air nominal'}];

hold on;
% plot([1:480] -0.5, PDD_vacuum_MC27./max(PDD_vacuum_MC27), '.-');
% leg = [leg, {'Vacuum MCemittance'}];

plot([1:480] -0.5 +2.05, PDD_air_MC27./max(PDD_air_MC27), '.-');
leg = [leg, {'Air MCemittance'}];

plot(machine.data(27).depths, PDD_basedata./max(PDD_basedata), '.-');
leg = [leg, {'machine'}];

grid on;
xlim([0, 110]);
legend(leg);
%legend('Vacuum nominal', 'Air nominal', 'Vacuum MCemittance', 'Air MCemittance', 'machine');
%% Same for carbon
integr = [80:80];
PDD_air_nominal31 = 'D:\matRad_gitHubRemo_MKM_TOPAS\PDD_carbon_air_nominal31.dcm';
PDD_air_nominal31 = squeeze(double(dicomread(PDD_air_nominal31)));
PDD_air_nominal31 = permute(PDD_air_nominal31, [2,3,1]);
PDD_air_nominal31 = squeeze(sum(PDD_air_nominal31(integr,integr,:), [1,2]));

PDD_vacuum_MC31 = 'D:\matRad_gitHubRemo_MKM_TOPAS\PDD_carbon_vacuum_MCemittance31.dcm';
PDD_vacuum_MC31 = squeeze(double(dicomread(PDD_vacuum_MC31)));
PDD_vacuum_MC31 = permute(PDD_vacuum_MC31, [2,3,1]);
PDD_vacuum_MC31 = squeeze(sum(PDD_vacuum_MC31(integr,integr,:), [1,2]));

PDD_air_MC31 = 'D:\matRad_gitHubRemo_MKM_TOPAS\PDD_carbon_air_MCemittance31.dcm';
PDD_air_MC31 = squeeze(double(dicomread(PDD_air_MC31)));
PDD_air_MC31 = permute(PDD_air_MC31, [2,3,1]);
PDD_air_MC31 = squeeze(sum(PDD_air_MC31(integr,integr,:), [1,2]));

%% plot
load('carbon_Generic.mat');
PDD_basedata = machine.data(31).Z;


leg = [];
figure;

hold on;
plot([1:480] -0.5, PDD_air_nominal31./max(PDD_air_nominal31), '.-');
leg = [leg, {'Air nominal'}];


hold on;
% plot([1:480] -0.5, PDD_vacuum_MC31./max(PDD_vacuum_MC31), '.-');
% leg = [leg, {'Vacuum MCemittance'}];
% 
plot([1:480] -0.5 + 2.05, PDD_air_MC31./max(PDD_air_MC31), '.-');
leg = [leg, {'Air MCemittance'}];

plot(machine.data(31).depths, PDD_basedata./max(PDD_basedata), '.-');
leg = [leg, {'machine'}];



grid on;
xlim([0, 110]);

legend(leg);