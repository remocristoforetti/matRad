%% Script: Read back TOPAS and compare with Hong for the helium/carbon/oxygen SOBP plans (entrance slabs)

%% set matRad runtime configuration
clear;
matRad_rc;
matRad_cfg = MatRad_Config.instance();

ions = {'helium','carbon','oxygen'};
nIons = numel(ions);

doseHong  = cell(1,nIons);
doseTopas = cell(1,nIons);

for i = 1:nIons
    % Load the plan information saved by script_ions_SOBPslabs_writeTopasFiles.m
    mainDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','doseCalculation','MCengines',[ions{i} 'SOBPslabs_2']);
    load(fullfile(mainDir,'planInfo.mat'),'ct','cst','pln','stf','w');

    % analytical (Hong) dose the weights were optimized with, as reference
    load(fullfile(matRad_cfg.primaryUserFolder,'dataPaper','doseCalculation','analyticalEngines',[ions{i} '_SOBPslabs_analytical_2.mat']),'doseSOBP');
    doseHong{i} = doseSOBP;

    % Read back TOPAS
    topasDir = fullfile(mainDir,'TOPAS');
    try
        engineTopas = DoseEngines.matRad_TopasMCEngine(pln);
        engineTopas.externalCalculation = topasDir;
        resultTopas  = engineTopas.calcDoseForward(ct,cst,stf,w);
        doseTopas{i} = resultTopas.physicalDose;
    catch ME
        matRad_cfg.dispWarning('%s: could not read TOPAS results from %s (%s)',ions{i},topasDir,ME.message);
    end

    % Save the results for later visualization
    doseHongIon  = doseHong{i};
    doseTopasIon = doseTopas{i};
    save(fullfile(matRad_cfg.primaryUserFolder,'dataPaper','doseCalculation','MCengines',[ions{i} '_SOBPslabs_MC_2.mat']), ...
        'ct','cst','pln','stf','w','doseHongIon','doseTopasIon');
end

%% Plot the depth-dose curves along the central beam axis (averaged over 3x3 central voxels)
isoIdx = matRad_world2cubeIndex(pln.propStf.isoCenter,ct); % [yIdx xIdx zIdx]
yIdx = isoIdx(1) + [-1 0 1];
zIdx = isoIdx(3) + [-1 0 1];

figure('Name','Hong vs TOPAS: ion SOBPs with entrance slabs','Color','w');
tiledlayout(1,nIons);
for i = 1:nIons
    nexttile;
    hold on;
    plot(ct.x,squeeze(mean(doseHong{i}(yIdx,:,zIdx),[1 3])),'LineWidth',2,'DisplayName','HongPB (optimized)');
    if ~isempty(doseTopas{i})
        plot(ct.x,squeeze(mean(doseTopas{i}(yIdx,:,zIdx),[1 3])),'LineWidth',2,'DisplayName','TOPAS');
    end
    grid on; box on;
    xlabel('x [mm]');
    ylabel('physical dose [Gy]');
    title(ions{i});
    legend('show','Location','northwest');
end

%% Plot the 2D dose distributions: Hong (top row) and TOPAS (bottom row) per ion
figure('Name','Hong vs TOPAS: ion SOBP dose distributions','Color','w');
tiledlayout(2,nIons);
for i = 1:nIons
    doseWindow = [0 max([doseHong{i}(:); doseTopas{i}(:)])];
    nexttile(i);
    matRad_plotSlice(ct,'axesHandle',gca,'cst',cst,'dose',doseHong{i}, ...
        'plane',3,'slice',isoIdx(3),'doseWindow',doseWindow,'title',[ions{i} ' - HongPB']);
    if ~isempty(doseTopas{i})
        nexttile(nIons + i);
        matRad_plotSlice(ct,'axesHandle',gca,'cst',cst,'dose',doseTopas{i}, ...
            'plane',3,'slice',isoIdx(3),'doseWindow',doseWindow,'title',[ions{i} ' - TOPAS']);
    end
end
