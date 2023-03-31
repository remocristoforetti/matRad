matRad_rc;
matRad_cfg = MatRad_Config.instance();
%% Setup dummy Phantom, %this is needed for how the topas config is setup but does not affact the scoring or simulation
resolution.x = 1; %mm
resolution.y = 1; %mm

resolution.z = 1; %mm
cubeRealDim = [100 100 100]; %mm
cubeDim  = ceil(cubeRealDim./[resolution.y resolution.x resolution.z]); %In Voxels

ct.cube = {ones(cubeDim)};
ct.resolution = resolution;
ct.cubeDim = cubeDim;
ct.numOfCtScen = 1;


ct.cubeHU = {zeros(cubeDim)}; %All water

ct.hlut = [1 0; 0 -1024 ];
ct.ctGrid.resolution = ct.resolution;

%% Setup

%pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = 0;

pln.propStf.couchAngles   = 0;
pln.propStf.bixelWidth    = 3;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
%pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

pln.radiationMode   = 'carbon';
pln.propDoseCalc.calcLET = 1;


%% Load machine
%Retrive some infos from Generic machine
%Generic_machine = load('carbon_Generic.mat');
HITmachine = load('carbon_HITgantry');
%machine = Generic_machine.machine;
machine = HITmachine.machine;
%% Dummy STF

%For the time being, keep the same energies saved in the Generic machine
%E = [machine.data(31:37).energy];
%E = [machine.data([31:10:64]).energy];
%E = [machine.data([78:10:112]).energy];
%E = [machine.data(31).energy];
%E = [machine.data(78).energy];
%E = [machine.data([78:2:142]).energy];
%E = [machine.data([31:64]).energy];
%E = [machine.data([31:40]).energy];
%E = [machine.data(:).energy];

E = [machine.data([1:20:end]).energy];
% Manually substitute 

% E = (E(1) + [-5:5])*12; %10 MeV range



stf.gantryAngle   = pln.propStf.gantryAngles;
stf.couchAngle    = pln.propStf.couchAngles;
stf.bixelWidth    = pln.propStf.bixelWidth;
stf.radiationMode = pln.radiationMode;
stf.SAD           = machine.meta.SAD;

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

%% Define structure for phantom info

[~, eIdx] = intersect([machine.data.energy], E);
phantom = cell(1,5); %cell(size(E,2),5);
for k=1:1%size(E,2)
    %Define eIdx
    phantom(k,1) = {eIdx(k)}; 

    %Define energy
    phantom(k,2) = {E(k)}; 

    %Define meta info
    phantom(k,3) = {struct('meta', [])};
    phantom{k,3}.meta.machine = 'generic';

    %Define energy binning parameters
    for ionIdx=1:size(Ions,2)
        phantom{k,4}(ionIdx) = struct( 'ion', [], ...
                                       'EMax', [], ...
                                       'EMin', [], ...
                                        'nEBins', []);
    end
    
    phantom{k,4}(1).ion = 'protons';
    phantom{k,4}(1).EMax = 500;
    phantom{k,4}(1).EMin = 0;
    phantom{k,4}(1).nEBins = 1000;

    phantom{k,4}(2).ion = 'He';
    phantom{k,4}(2).EMax = 400;
    phantom{k,4}(2).EMin = 0;
    phantom{k,4}(2).nEBins = 800;

    phantom{k,4}(3).ion = 'Li';
    phantom{k,4}(3).EMax = 400;
    phantom{k,4}(3).EMin = 0;
    phantom{k,4}(3).nEBins = 800;

    phantom{k,4}(4).ion = 'Be';
    phantom{k,4}(4).EMax = 400;
    phantom{k,4}(4).EMin = 0;
    phantom{k,4}(4).nEBins = 800;

    phantom{k,4}(5).ion = 'B';
    phantom{k,4}(5).EMax = 300;
    phantom{k,4}(5).EMin = 0;
    phantom{k,4}(5).nEBins = 600;

    phantom{k,4}(6).ion = 'C';
    phantom{k,4}(6).EMax = 300;
    phantom{k,4}(6).EMin = 0;
    phantom{k,4}(6).nEBins = 600;

    %Define phantom specifications
    phantom{k,5}.resolution.x = 5;
    phantom{k,5}.resolution.y = 0.1;
    phantom{k,5}.resolution.z = 5;
    
    phantom{k,5}.name = 'Phantom';
    phantom{k,5}.dimension.x = 100; %in mm
    phantom{k,5}.dimension.y = 350; %in mm
    phantom{k,5}.dimension.z = 100; %in mm
    


    phantom{k,5}.includePDD = true;
    phantom{k,5}.includeDS = false;
    phantom{k,5}.includeSP = false;

    phantom{k,5}.includeEdEventSpectrum = false;
    phantom{k,5}.includeProtonLET = false;
    phantom{k,5}.includeSPdose = false;
    phantom{k,5}.includeEventsCounter = false;

end





for k=1:size(phantom{1,4},2)
 phantom{1,4}(k).EBinWidth = (phantom{1,4}(k).EMax-phantom{1,4}(k).EMin)/phantom{1,4}(k).nEBins;
 phantom{1,4}(k).E = linspace(phantom{1,4}(k).EMin,...
                                                   phantom{1,4}(k).EMax - phantom{1,4}(k).EBinWidth,...
                                                   phantom{1,4}(k).nEBins)...
                                                   +phantom{1,4}(k).EBinWidth/2;
end
%% Setup topas config

topas_cfg = matRad_TopasConfig;
topas_cfg.worldMaterial = 'Vacuum';

folderName = 'BaseData_carbon_spectra\PDDsEsFirst';

if ~exist(strcat(topas_cfg.workingDir, folderName), 'dir')
   mkdir(strcat(topas_cfg.workingDir, folderName));
end

topas_cfg.workingDir = strcat(topas_cfg.workingDir, folderName);
topas_cfg.label = 'BaseData';
topas_cfg.numOfRuns = 1;

topas_cfg.numThreads = -1;

numHistoriesPerBixel = 1000;

topas_cfg.numHistories = ceil(numHistoriesPerBixel*size(stf.ray.energy,2));
topas_cfg.scorer.outputType      = 'DICOM';
topas_cfg.scorer.reportQuantity  = {'Sum'};
topas_cfg.scorer.volume          = true;
topas_cfg.scorer.doseToMedium    = false; %Turn down and add parallel scorers for higher resolution

topas_cfg.scorer.filename        = 'constructor';
w = (1/size(stf.ray.energy,2))*ones(size(stf.ray.energy,2),1);

writeScoringTreeDirectory_singlePhantom(phantom, Ions,topas_cfg);
topas_cfg.writeAllFiles(ct,0,pln,stf,machine,w);
%% Data analisis
analisis = highResPhantomAnalisis_singlePhantom();

%analisis.wDir     = 'C:\r408i_data\r408i_data\BaseData_carbon_spectra\NewScorers'; %topas_cfg.workingDir;
%analisis.wDir     = 'C:\r408i_data\r408i_data\BaseData_carbon_spectra\NewScorers\HITgantry\SOBP';
analisis.wDir      = 'C:\r408i_data\r408i_data\BaseData_carbon_spectra\NewScorers\SOBP_Generic_cluster';

analisis.readPDD   = true;

analisis.readSpectra = true; % This is fluence scorer
analisis.readDSSpectra = false;
analisis.readSPeDSpectra = true;
analisis.readProtonLET = false;

analisis.readSPdose = true;
analisis.readEventsCounter = true;

analisis.ions = {'protons', 'He', 'Li', 'Be', 'B', 'C'};
analisis.phantom  = phantom;
analisis.integrationSurface = [1:20];
%analisis.averageQuantities = false;

analisis.averageQuantities = true; % true if using counter scorers

analisis.downSamplingPhantom.resolutions = [1 0.1 3]; %[1 0.1 3];
analisis.downSamplingPhantom.highResWindowWidth = 30;  %in mm

%Set eLine(2) to -1 for automatic sampling
% analisis.downSamplingPhantom.EParam(1).nEnergies = 1000;
% analisis.downSamplingPhantom.EParam(1).eLine = [phantom{1,4}(1).E(1) 0.499]; % eMax = eLine(2)*nEnergies + eLine(1); constant for all depths, then modify and make id depth dependent -> for every downsampled depth
% 
% analisis.downSamplingPhantom.EParam(2).nEnergies = 800;
% 
% analisis.downSamplingPhantom.EParam(2).eLine = [phantom{1,4}(2).E(1) 0.499]; % eMax = eLine(2)*nEnergies + eLine(1); constant for all depths, then modify and make id depth dependent -> for every downsampled depth
% 
% analisis.downSamplingPhantom.EParam(3).nEnergies = 800;
% analisis.downSamplingPhantom.EParam(3).eLine = [phantom{1,4}(3).E(1) 0.499]; % eMax = eLine(2)*nEnergies + eLine(1); constant for all depths, then modify and make id depth dependent -> for every downsampled depth
% 
% analisis.downSamplingPhantom.EParam(4).nEnergies = 800;
% analisis.downSamplingPhantom.EParam(4).eLine = [phantom{1,4}(4).E(1) 0.499]; % eMax = eLine(2)*nEnergies + eLine(1); constant for all depths, then modify and make id depth dependent -> for every downsampled depth
% 
% analisis.downSamplingPhantom.EParam(5).nEnergies = 600;
% 
% analisis.downSamplingPhantom.EParam(5).eLine = [phantom{1,4}(5).E(1) 0.499]; % eMax = eLine(2)*nEnergies + eLine(1); constant for all depths, then modify and make id depth dependent -> for every downsampled depth
% 
% analisis.downSamplingPhantom.EParam(6).nEnergies = 1000;
% analisis.downSamplingPhantom.EParam(6).eLine = [phantom{1,4}(6).E(1) 0.299]; % eMax = eLine(2)*nEnergies + eLine(1); constant for all depths, then modify and make id depth dependent -> for every downsampled depth

%Set eLine(2) to -1 for automatic sampling

analisis.downSamplingPhantom.EParam(1).nEnergies = 250;
analisis.downSamplingPhantom.EParam(1).eLine = [phantom{1,4}(1).E(1) 1]; % eMax = eLine(2)*nEnergies + eLine(1); constant for all depths, then modify and make id depth dependent -> for every downsampled depth

analisis.downSamplingPhantom.EParam(2).nEnergies = 250;
analisis.downSamplingPhantom.EParam(2).eLine = [phantom{1,4}(2).E(1) 1]; % eMax = eLine(2)*nEnergies + eLine(1); constant for all depths, then modify and make id depth dependent -> for every downsampled depth

analisis.downSamplingPhantom.EParam(3).nEnergies = 250;
analisis.downSamplingPhantom.EParam(3).eLine = [phantom{1,4}(3).E(1) 1]; % eMax = eLine(2)*nEnergies + eLine(1); constant for all depths, then modify and make id depth dependent -> for every downsampled depth

analisis.downSamplingPhantom.EParam(4).nEnergies = 250;
analisis.downSamplingPhantom.EParam(4).eLine = [phantom{1,4}(4).E(1) 1]; % eMax = eLine(2)*nEnergies + eLine(1); constant for all depths, then modify and make id depth dependent -> for every downsampled depth

analisis.downSamplingPhantom.EParam(5).nEnergies = 250;
analisis.downSamplingPhantom.EParam(5).eLine = [phantom{1,4}(5).E(1) 1]; % eMax = eLine(2)*nEnergies + eLine(1); constant for all depths, then modify and make id depth dependent -> for every downsampled depth

analisis.downSamplingPhantom.EParam(6).nEnergies = 500;
analisis.downSamplingPhantom.EParam(6).eLine = [phantom{1,4}(6).E(1) 0.499]; % eMax 



bD = [];
for k=[1,3,26,32]%k=1:size(E,2)
    analisis.run       = k-1;
    analisis.importRawData();


    DS(k).Phi       = analisis.Phi;
    DS(k).PDD       = analisis.PDD;
    DS(k).depths    = analisis.depths;
    DS(k).edPhi     = analisis.edPhi;
    DS(k).fPhi      = analisis.fPhi;

    %    DS(k).ProtonLET = analisis.ProtonLET;

    DS(k).newPhi       = analisis.newPhi;
    DS(k).newPDD       = analisis.newPDD;
    
    DS(k).newDepths    = analisis.newDepths;
    DS(k).newedPhi     = analisis.newedPhi;
    DS(k).eBinning  = analisis.EParam(:);
    DS(k).spectraLET = analisis.spectraLET;
    DS(k).spectraLETdS = analisis.spectraLETdS; % downsampled
    DS(k).newfPhi =     analisis.newfPhi;
    DS(k).spectraLETt = analisis.spectraLETt;
    DS(k).dosePhi = analisis.dosePhi;
    
    for m=1:6
        %DS(k).EParam(m).E = (DS(k).eBinning(m).eLine(:,2).*[0:DS(k).eBinning(m).nEnergies-1] + DS(k).eBinning(m).eLine(:,1))';
        DS(k).E{m} = phantom{1,4}(m).E;
    end
    matRad_cfg.dispInfo('analizing E %u \n',k);

    bD(k).PDD       = analisis.newPDD;
%    bD(k).Phi       = analisis.newPhi;
    bD(k).fPhi       = analisis.newfPhi;

    
    bD(k).edPhi     = analisis.newedPhi;

    
    bD(k).depths    = analisis.newDepths;
    bD(k).eBinning  = analisis.EParam(:);
%    bD(k).dosePhi   = analisis.newDosePhi;

    bD(k).PDD       = analisis.PDD; 


    bD(k).Phi       = analisis.Phi;
    bD(k).fPhi       = analisis.fPhi;

    bD(k).edPhi     = analisis.edPhi;


    bD(k).depths    = analisis.depths;
    bD(k).eBinning = DS(k).E;
%      for m=1:6

%          bD(k).eBinning(m) = {phantom{1,4}(m).E};
%      end
end

%% NaN

for k=1:size(eIdx,1)
    for m=1:6
        
        v = sum(isnan(bD(k).Phi{m}), [1,2,3]);

        if v>0
            display(['NaN: k = ', num2str(k), ' m = ', num2str(m)]);
        end
    end

    display('done');
end

%% Load into Base Data

newMachine = machine;

[~,eIdx] = intersect([machine.data.energy], E);


for k=[1,3,26,32]%k=1:size(eIdx,1)

    newMachine.data(eIdx(k)).DS = bD(k);
end

% 
% oldMachine = machine;
% machine = newMachine;
% save(['C:/r408i_data/r408i_data/carbon_newGenericSOBP_FullRes.mat'], 'machine');
% machine = oldMachine;

%% PDD comparison
[~,eIdx] = intersect([machine.data.energy], E);
selectedEnergy = 21;
figure;
plot(bD(selectedEnergy).depths, bD(selectedEnergy).PDD./max(bD(selectedEnergy).PDD), '.-');
hold on;
plot(machine.data(eIdx(selectedEnergy)).depths, machine.data(eIdx(selectedEnergy)).Z./max(machine.data(eIdx(selectedEnergy)).Z), '.-');
grid on;
grid minor;

legend('Simulation', 'BaseData', 'FontSize', 14);
xlabel('Depth [mm]', 'FontSize',14);
ylabel('Dose', 'FontSize', 14);
%% Plot some info

ionSelect = 6;
eIdx = 1;
Phi = full(DS(eIdx).Phi{ionSelect});
newPhi = full(DS(eIdx).newPhi{ionSelect});

edPhi = full(DS(eIdx).edPhi{ionSelect});

newedPhi = full(DS(eIdx).newedPhi{ionSelect});

depthSelect = 95; %in mm

depthSelectIdx = interp1(DS(eIdx).depths,[1:size(DS(eIdx).depths,2)],depthSelect, 'nearest');
newDepthSelectIdx = interp1(DS(eIdx).newDepths,[1:size(DS(eIdx).newDepths,2)], depthSelect, 'nearest');
prof= Phi(:,depthSelectIdx);
newProf = newPhi(:,newDepthSelectIdx);

profeD= edPhi(:,depthSelectIdx);
profeD_s1= edPhi(:,depthSelectIdx);

newProfeD = newedPhi(:,newDepthSelectIdx);
EBinWidth = (phantom{1,4}(ionSelect).EMax-phantom{1,4}(ionSelect).EMin)/phantom{1,4}(ionSelect).nEBins;
eBin = linspace(phantom{1,4}(k).EMin,...
               phantom{1,4}(ionSelect).EMax - EBinWidth,...
               phantom{1,4}(ionSelect).nEBins)...
               +EBinWidth/2;

newEBin = DS(eIdx).EParam(ionSelect).E(:,newDepthSelectIdx);

% figure;
% plot(eBin, prof, '.-');
% hold on;
% plot(newEBin,newProf, '.-');
% legend('old', 'new');
% grid on;
% grid minor;

figure;
plot(eBin, profeD, '.-');
hold on;
plot(newEBin,newProfeD, '.-');
legend('old', 'new');

grid on;
grid minor;
%% plot comparison Phi fPhi
ionSelect = 6;
eIdx = 1;

Phi = full(DS(eIdx).Phi{ionSelect});
fPhi = full(DS(eIdx).fPhi{ionSelect});

depthSelect = 95.7; %in mm

depthSelectIdx = interp1(DS(eIdx).depths,[1:size(DS(eIdx).depths,2)],depthSelect, 'nearest');
fDepthSelectIdx = interp1(DS(eIdx).depths,[1:size(DS(eIdx).depths,2)], depthSelect, 'nearest');
prof= Phi(:,depthSelectIdx);
fProf = fPhi(:,fDepthSelectIdx);

EBinWidth = (phantom{1,4}(ionSelect).EMax-phantom{1,4}(ionSelect).EMin)/phantom{1,4}(ionSelect).nEBins;

eBin = linspace(phantom{1,4}(k).EMin,...
               phantom{1,4}(ionSelect).EMax - EBinWidth,...
               phantom{1,4}(ionSelect).nEBins)...
               +EBinWidth/2;

%newEBin = DS(eIdx).EParam(ionSelect).E(:,newDepthSelectIdx);
figure;

plot(eBin, prof./max(prof), '.-');
hold on;

plot(eBin,fProf./max(fProf), '.-');
%plot(newEBin,newProfeD./max(newProfeD), '.-');

legend('Phi', 'fPhi');
grid on;
grid minor;


%% Add BB let calculation for comparison with baseData
[~,eIdx] = intersect([machine.data.energy], E);

selectedEnergy = 1;
figure;
plot(machine.data(eIdx(selectedEnergy)).depths,machine.data(eIdx(selectedEnergy)).LET, '.-');
hold on;

spectraLET =  zeros(size(DS(selectedEnergy).spectraLET{1}));
spectraLETt = zeros(size(DS(selectedEnergy).spectraLET{1}));


for k=1:1
    spectraLET = spectraLET + DS(selectedEnergy).spectraLET{k};
    spectraLETt = spectraLETt + DS(selectedEnergy).spectraLETt{k};
end
plot(DS(selectedEnergy).depths, spectraLET, '.-');
% hold on;
% plot(DS(selectedEnergy).depths, spectraLETt, '.-');

legend('base Data', 'simulated', 'FontSize', 14);
grid on;

grid minor;
xlabel('Depths [mm]', 'FontSize', 14);
ylabel('LET keV/mu', 'FontSize', 14);

%% Compare PDD and Spectra

k=1;

figure;
movegui(['northwest']);
plot(machine.data(eIdx(k)).depths, machine.data(eIdx(k)).Z./max(machine.data(eIdx(k)).Z), '.-');
hold on;
plot(DS(k).depths, DS(k).PDD./max(DS(k).PDD), '.-');

