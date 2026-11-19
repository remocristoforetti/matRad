%% Script: Read back and compare TOPAS/MCsquare/FRED on the optimized SOBP (entrance slabs) plan

%% set matRad runtime configuration
clear;
matRad_rc;

matRad_cfg = MatRad_Config.instance();

%% Load the plan information saved by script_protons_SOBPslabs_writeParameterFiles.m
mainDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','doseCalculation','MCengines','protonsSOBPslabs_2');
load(fullfile(mainDir,'planInfo.mat'),'ct','cst','pln','stf','w');

% analytical (Hong) dose the weights were optimized with, as reference
load(fullfile(matRad_cfg.primaryUserFolder,'dataPaper','doseCalculation','analyticalEngines','protons_SOBPslabs_analytical_2.mat'),'doseSOBP');
doseHong = doseSOBP;

%% Read back the three MC simulations
topasDir    = fullfile(mainDir,'TOPAS');    % folder holding MCparam.mat and the score_*.bin files
mcSquareDir = fullfile(mainDir,'MCsquare'); % engine reads <folder>/output
fredDir     = fullfile(mainDir,'FRED');     % engine appends MCrun/out itself

doseTopas    = [];
doseMCsquare = [];
doseFred     = [];

% TOPAS
engineTopas = DoseEngines.matRad_TopasMCEngine(pln);
try
    engineTopas.externalCalculation = topasDir;
    resultTopas = engineTopas.calcDoseForward(ct,cst,stf,w);
    doseTopas   = resultTopas.physicalDose;
catch ME
    matRad_cfg.dispWarning('TOPAS: could not read results from %s (%s)',topasDir,ME.message);
end

% MCsquare
engineMCsquare = DoseEngines.matRad_ParticleMCsquareEngine(pln);
try
    engineMCsquare.workingDir          = mcSquareDir; % keep internal file paths consistent with the read folder
    engineMCsquare.externalCalculation = mcSquareDir;
    resultMCsquare = engineMCsquare.calcDoseForward(ct,cst,stf,w);
    doseMCsquare   = resultMCsquare.physicalDose;
catch ME
    matRad_cfg.dispWarning('MCsquare: could not read results from %s (%s)',fullfile(mcSquareDir,'output'),ME.message);
end

% FRED
engineFred = DoseEngines.matRad_ParticleFREDEngine(pln);
try
    engineFred.externalCalculation = fredDir;
    resultFred = engineFred.calcDoseForward(ct,cst,stf,w);
    doseFred   = resultFred.physicalDose;
catch ME
    matRad_cfg.dispWarning('FRED: could not read results from %s (%s)',fullfile(fredDir,'MCrun','out'),ME.message);
end

engineNames = {'HongPB (optimized)', engineTopas.name, engineMCsquare.name, engineFred.name};

%% Save the results for later visualization
outFile = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','doseCalculation','MCengines','protons_SOBPslabs_MC_2.mat');
save(outFile,'ct','cst','pln','stf','w', ...
    'doseHong','doseTopas','doseMCsquare','doseFred','engineNames');
matRad_cfg.dispInfo('Results saved to %s\n',outFile);

%% Keep only engines with actual (non-empty, non-zero) dose for plotting

allDoses = {doseHong, doseTopas, doseMCsquare, doseFred};
hasDose  = cellfun(@(d) ~isempty(d) && max(d(:)) > 0,allDoses);
doses    = allDoses(hasDose);
names    = engineNames(hasDose);
nEngines = numel(doses);
matRad_cfg.dispInfo('Plotting %d of %d engines: %s\n',nEngines,numel(allDoses),strjoin(names,', '));

%% Plot the depth-dose curves along the central beam axis
isoIdx = matRad_world2cubeIndex(pln.propStf.isoCenter,ct); % [yIdx xIdx zIdx]
centralProfile(1,:) = isoIdx(1) + [-1,0,1];
centralProfile(2,:) = isoIdx(3) + [-1,0,1];

figure('Name','MC engine comparison: SOBP slabs','Color','w');
hold on;
for i = 1:nEngines
    profile = squeeze(mean(doses{i}(centralProfile(1,:),:,centralProfile(2,:)), [1,3]));
    plot(ct.x,profile,'LineWidth',2,'DisplayName',names{i});
end
xline(pln.propStf.isoCenter(1),'k--','DisplayName','Target center');
grid on; box on;
xlabel('x [mm]');
ylabel('physical dose [Gy]');
title('Optimized proton SOBP with entrance slabs: depth dose along the central beam axis');
legend('show','Location','northwest');

%% Plot the 2D dose distributions side by side for comparison
doseWindow = [0 max(cellfun(@(d) max(d(:)),doses))];

figure('Name','MC engine comparison: SOBP slabs dose distributions','Color','w');
for i = 1:nEngines
    nexttile;
    matRad_plotSlice(ct,'axesHandle',gca,'cst',cst,'dose',doses{i}, ...
        'plane',3,'slice',isoIdx(3),'doseWindow',doseWindow,'title',names{i});
end
