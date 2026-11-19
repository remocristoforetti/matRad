%% Script 2: Compare analytical proton pencil-beam engines on a single PB
%
% Loads the box phantom built in
% scriptsForPaper/phantomBuilder/matRad_script01_buildPhantom.m (saved as
% phantom_insert.mat) and fires a single proton pencil beam through it,
% once per analytical dose-calculation engine (Hong, FineSampling,
% AnalyticalBortfeld). Each engine is instantiated by hand as its own
% object - bypassing the automatic pln.propDoseCalc.engine selection used
% by matRad_calcDoseForward - so that engine-specific options can be set
% explicitly. The resulting central-axis depth-dose curves are then
% compared in a single plot.
%
% NOTE: phantom_insert.mat contains an extra 'Insert' OAR between the
% beam entrance and the target, presumably meant as a density
% heterogeneity. However, matRad_PhantomBuilder.getctcst() assigns HU in
% reverse VOI-definition order (first-added VOI wins overlaps), and
% 'PHANTOM' - added first and spanning the whole cube - overwrites
% 'Insert' everywhere they overlap. The loaded ct is therefore uniform
% water (HU = 0) despite the extra structure; this script compares the
% three engines' baseline behavior in homogeneous water, not their
% heterogeneity handling.

%% set matRad runtime configuration
clear;
matRad_rc;
matRad_cfg = MatRad_Config.instance();

%% Load the phantom built in phantomBuilder/matRad_script01_buildPhantom.m
load(fullfile(matRad_cfg.primaryUserFolder,'dataPaper','phantomBuilder','phantom_insert.mat'),'ct','cst');

%% Treatment plan meta information
pln.radiationMode = 'protons';
pln.machine       = 'Generic';
pln.bioModel      = 'none';    % pure physical dose comparison
pln.multScen      = 'nomScen'; % single nominal scenario, no robustness

% Beam entering from the -x face and travelling in +x, towards the offset
% target. (matRad_getRotationMatrix: at gantry 270 deg the source sits at
% -x and the beam travels towards +x - see script 1 for the derivation.)
pln.propStf.gantryAngles = 270;
pln.propStf.couchAngles  = 0;
pln.propStf.bixelWidth   = 5; % [mm], irrelevant for a single beamlet, kept for a valid stf generator setup
pln.propStf.isoCenter    = matRad_getIsoCenter(cst,ct,0); % target centroid

% match the dose grid resolution to the (1 mm) ct resolution for a smooth
% depth-dose curve
pln.propDoseCalc.doseGrid.resolution.x = 1;
pln.propDoseCalc.doseGrid.resolution.y = 1;
pln.propDoseCalc.doseGrid.resolution.z = 1;

%% Generate steering information for a single pencil beam
% matRad_StfGeneratorParticleSingleBeamlet places exactly one spot, aimed
% at pln.propStf.isoCenter, and automatically picks the energy whose
% Bragg peak lands there (ray-traced radiological depth to the isocenter).
stfGenerator = matRad_StfGeneratorParticleSingleBeamlet(pln);
stf = stfGenerator.generate(ct,cst);

w = 1; % single bixel, unit weight -> "one PB"

%% Instantiate the three analytical dose engines explicitly
% Each engine is its own object built from the same pln, so
% engine-specific properties can be inspected/overridden by hand. This is
% what DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln) would do
% for you automatically inside matRad_calcDoseForward/matRad_calcDoseInfluence.
engineHong = DoseEngines.matRad_ParticleHongPencilBeamEngine(pln);
engineHong.lateralModel = 'single'; % force a single-Gaussian lateral optics model (default: 'fast')

engineFineSampling = DoseEngines.matRad_ParticleFineSamplingPencilBeamEngine(pln);
engineFineSampling.fineSampling.N = 2; % more sub-beam samples per spot (default: 2) -> more accurate lateral sampling

engineBortfeld = DoseEngines.matRad_ParticleAnalyticalBortfeldEngine(pln);
engineBortfeld.sigmaEnergy = 0.02; % widen the assumed Gaussian energy spread of the beam (default: 0.01) -> softer distal fall-off

%% Forward dose calculation with each engine
resultHong          = engineHong.calcDoseForward(ct,cst,stf,w);
resultFineSampling  = engineFineSampling.calcDoseForward(ct,cst,stf,w);
resultBortfeld      = engineBortfeld.calcDoseForward(ct,cst,stf,w);

%% Save the results for later visualization
% Only the physical dose cubes are kept (not the full result structs),
% together with the engine names for plot labels. ct and pln are saved
% as well since plotting needs the ct grid and the isocenter.
doseHong          = resultHong.physicalDose;
doseFineSampling  = resultFineSampling.physicalDose;
doseBortfeld      = resultBortfeld.physicalDose;
engineNames       = {engineHong.name, engineFineSampling.name, engineBortfeld.name};

outDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','doseCalculation','analyticalEngines');
if ~isfolder(outDir)
    mkdir(outDir);
end
outFile = fullfile(outDir,'protons_insert_analytical.mat');
save(outFile,'ct','cst','pln','stf','w', ...
    'doseHong','doseFineSampling','doseBortfeld','engineNames');
%%

outDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','doseCalculation','analyticalEngines');
outFile = fullfile(outDir,'protons_insert_analytical.mat');

matRad_cfg = MatRad_Config.instance();
load(fullfile(outFile));
%% Extract the central-axis depth-dose profile for each engine
% The beam travels along x through the isocenter's y/z position. Since
% the dose grid resolution was set to match the ct grid, ct-based indices
% can be used directly to index into the dose cubes.
isoIdx = matRad_world2cubeIndex(pln.propStf.isoCenter,ct); % [yIdx xIdx zIdx]

profileHong         = squeeze(doseHong(isoIdx(1),:,isoIdx(3)));
profileFineSampling = squeeze(doseFineSampling(isoIdx(1),:,isoIdx(3)));
profileBortfeld     = squeeze(doseBortfeld(isoIdx(1),:,isoIdx(3)));

% %% Plot the three depth-dose curves for comparison
% figure('Name','Analytical engine comparison: single PB','Color','w');
% hold on;
% plot(ct.x,profileHong,'LineWidth',2,'DisplayName',engineHong.name);
% plot(ct.x,profileFineSampling,'LineWidth',2,'DisplayName',engineFineSampling.name);
% plot(ct.x,profileBortfeld,'LineWidth',2,'DisplayName',engineBortfeld.name);
% xline(pln.propStf.isoCenter(1),'k--','DisplayName','Target center');
% grid on; box on;
% xlabel('x [mm]');
% ylabel('physical dose [Gy] (per unit spot weight)');
% title('Single proton pencil beam: depth dose along the central beam axis');
% legend('show','Location','northwest');

%% Plot the three 2D dose distributions side by side for comparison
% Same axial slice (through the isocenter) and the same dose window for
% all three engines, so the color scale is directly comparable.
doseWindow = [0,1.8e-3];%[0 max([resultHong.physicalDose(:); resultFineSampling.physicalDose(:); resultBortfeld.physicalDose(:)])];

figure('Name','Analytical engine comparison: dose distributions','Color','w');

nexttile;
matRad_plotSlice(ct,'axesHandle',gca,'cst',cst,'dose',doseHong, ...
    'plane',3,'slice',isoIdx(3),'doseWindow',doseWindow,'title','Hong PB');

nexttile;
matRad_plotSlice(ct,'axesHandle',gca,'cst',cst,'dose',doseFineSampling, ...
    'plane',3,'slice',isoIdx(3),'doseWindow',doseWindow,'title','Fine Sampling');

nexttile;

matRad_plotSlice(ct,'axesHandle',gca,'cst',cst,'dose',doseBortfeld, ...
    'plane',3,'slice',isoIdx(3),'doseWindow',doseWindow,'title','Bortfeld');
%% PLOT TO TICKZ
% saveDir = fullfile(matRad_cfg.primaryUserFolder, 'dataPaper', 'Figures', 'doseCalculation', 'analyticalEngines', 'protons_insert', 'Hong');
% saveDir = fullfile(matRad_cfg.primaryUserFolder, 'dataPaper', 'Figures', 'doseCalculation', 'analyticalEngines', 'protons_insert', 'fineSampling');
saveDir = fullfile(matRad_cfg.primaryUserFolder, 'dataPaper', 'Figures', 'doseCalculation', 'analyticalEngines', 'protons_insert', 'Bortfeld');

cstPlot = cst;
% doseCube = doseHong;
% doseCube = doseFineSampling;
doseCube = doseBortfeld;

savePlots = false;

lowCTLim =-400;
higCtLim =1140;

level = -580;

f = figure;
movegui(gcf(), 'southwest');
f.Position(1:2) = f.Position(1:2);%-800;
f.Position(3:4) = [900,900];
[hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSlice(ct, 'cst', cstPlot, ...
                 'axesHandle', gca(),...
                  'dose', double(doseCube), ...
                  'plane', 3, ...
                  'slice',27,...
                  'doseIsoLevels', 1.8e-3*[50,70,80,95,105]./100,...
                  'ctWindow', [lowCTLim,higCtLim]+level, ...
                  'alpha', 0.55,...
                  'doseWindow', [0,1.8e-3],...
                  'voiSelection', [1 1 0]);...
                  % 'doseColorMap', diffMap);

% xlim([75,155]-5);
% ylim([75,155]-5);

title('');

if savePlots
    fileName = sprintf('DoseDist.tex');
    
    figSaveDir = fullfile(saveDir);
    
    if ~exist(figSaveDir, 'dir')
        mkdir(figSaveDir);
    end
    
    matRad_saveFigToTikz(f, figSaveDir, fileName, hContour, hIsoDose);
end