%% Script: Read back and compare TOPAS/MCsquare/FRED on a single PB

%% set matRad runtime configuration
clear;
matRad_rc;

matRad_cfg = MatRad_Config.instance();

%% Load the plan information saved by script_protons_insert_writeParameterFiles.m
mainDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','doseCalculation','MCengines','protonsInsert');
load(fullfile(mainDir,'planInfo.mat'),'ct','cst','pln','stf','w');

%% Read back the three MC simulations
topasDir    = fullfile(mainDir,'TOPAS');    % folder holding MCparam.mat and the score_*.bin files
mcSquareDir = fullfile(mainDir,'MCsquare'); % engine reads <folder>/output
fredDir     = fullfile(mainDir,'FRED');     % engine appends MCrun/out itself

doseTopas    = [];
doseMCsquare = [];
doseFred     = [];

% TOPAS
engineTopas = DoseEngines.matRad_TopasMCEngine(pln);
engineTopas.externalCalculation = topasDir;
resultTopas = engineTopas.calcDoseForward(ct,cst,stf,w);
doseTopas   = resultTopas.physicalDose;

% MCSquare
engineMCsquare = DoseEngines.matRad_ParticleMCsquareEngine(pln);
engineMCsquare.workingDir          = mcSquareDir; % keep internal file paths consistent with the read folder
engineMCsquare.externalCalculation = mcSquareDir;
resultMCsquare = engineMCsquare.calcDoseForward(ct,cst,stf,w);
doseMCsquare   = resultMCsquare.physicalDose;

% FRED
engineFred = DoseEngines.matRad_ParticleFREDEngine(pln);
engineFred.externalCalculation = fredDir;
resultFred = engineFred.calcDoseForward(ct,cst,stf,w);
doseFred   = resultFred.physicalDose;

engineNames = {engineTopas.name, engineMCsquare.name, engineFred.name};

%% Save the results for later visualization
% Same layout as the analytical-engine comparison. Engines that could not
% be read are saved as empty arrays.
outFile = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','doseCalculation','MCengines','protons_insert_MC.mat');
save(outFile,'ct','cst','pln','stf','w', ...
    'doseTopas','doseMCsquare','doseFred','engineNames');
matRad_cfg.dispInfo('Results saved to %s\n',outFile);

%% Keep only engines with actual (non-empty, non-zero) dose for plotting
allDoses = {doseTopas, doseMCsquare, doseFred};
hasDose  = cellfun(@(d) ~isempty(d) && max(d(:)) > 0,allDoses);
doses    = allDoses(hasDose);
names    = engineNames(hasDose);
nEngines = numel(doses);

if nEngines == 0
    matRad_cfg.dispError('No MC engine has results to plot yet - run script_protons_insert_runMCEngines.m first.');
end
matRad_cfg.dispInfo('Plotting %d of %d engines: %s\n',nEngines,numel(allDoses),strjoin(names,', '));

%% Plot the depth-dose curves along the central beam axis
% The beam travels along x through the isocenter's y/z position; the dose
% grid matches the ct grid, so ct-based indices index the dose cubes.
isoIdx = matRad_world2cubeIndex(pln.propStf.isoCenter,ct); % [yIdx xIdx zIdx]

figure('Name','MC engine comparison: single PB','Color','w');
hold on;
for i = 1:nEngines
    profile = squeeze(doses{i}(isoIdx(1),:,isoIdx(3)));
    plot(ct.x,profile,'LineWidth',2,'DisplayName',names{i});
end
xline(pln.propStf.isoCenter(1),'k--','DisplayName','Target center');
grid on; box on;
xlabel('x [mm]');
ylabel('physical dose [Gy] (per unit spot weight)');
title('Single proton pencil beam: depth dose along the central beam axis');
legend('show','Location','northwest');

%% Plot the 2D dose distributions side by side for comparison
% Same axial slice (through the isocenter) and the same dose window for
% all engines, so the color scale is directly comparable.
doseWindow = [0 max(cellfun(@(d) max(d(:)),doses))];

figure('Name','MC engine comparison: dose distributions','Color','w');
for i = 1:nEngines
    nexttile;
    matRad_plotSlice(ct,'axesHandle',gca,'cst',cst,'dose',doses{i}, ...
        'plane',3,'slice',isoIdx(3),'doseWindow',doseWindow,'title',names{i});
end
