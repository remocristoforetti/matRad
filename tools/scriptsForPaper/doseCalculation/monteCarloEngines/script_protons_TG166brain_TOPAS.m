%% Script: Proton plan on the TG166 brain phantom (HongPB) and TOPAS recalculation
% Optimizes a physical-dose proton plan with HongPB, then recalculates the
% optimized weights with TOPAS. Three TOPAS steps, each toggled below:
% write the input files, run TOPAS, read the results back into matRad.
clear;
matRad_rc;
matRad_cfg = MatRad_Config.instance();

optimizePlan = true; % false: load the previously optimized plan
writeFiles   = true;
runTopas     = true;
readTopas    = true;

numHistoriesDirect = 1e6;

topasExecCommand = 'wsl export TOPAS_G4_DATA_DIR=~/Applications/GEANT4/G4DATA; ~/Applications/TOPAS/OpenTOPAS-install/bin/topas';

mainDir  = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','doseCalculation','MCengines','protonsTG166brain');
topasDir = fullfile(mainDir,'TOPAS');
if ~isfolder(topasDir)
    mkdir(topasDir);
end
planFile = fullfile(mainDir,'planInfo.mat');

%% HongPB optimization
if optimizePlan
    load(fullfile(matRad_cfg.primaryUserFolder,'patients','TG166_brain_ct3mm.mat'),'ct','cst');

    % the phantom only ships with target objectives, set all objectives here
    cst(:,6) = {[]};
    getStructIdx = @(name) find(strcmp(cst(:,2),name));

    prescription = 60; % [Gy], total over all fractions

    cst{getStructIdx('PTV'),6}{1}          = DoseObjectives.matRad_SquaredDeviation(800,prescription);
    cst{getStructIdx('GTV'),6}{1}          = DoseObjectives.matRad_SquaredDeviation(800,prescription);
    cst{getStructIdx('Brain Stem'),6}{1}   = DoseObjectives.matRad_SquaredOverdosing(300,50);
    cst{getStructIdx('OpticChiasm'),6}{1}  = DoseObjectives.matRad_SquaredOverdosing(300,50);
    cst{getStructIdx('R OpticNerve'),6}{1} = DoseObjectives.matRad_SquaredOverdosing(300,50);
    cst{getStructIdx('L OpticNerve'),6}{1} = DoseObjectives.matRad_SquaredOverdosing(300,50);
    cst{getStructIdx('R Eye'),6}{1}        = DoseObjectives.matRad_SquaredOverdosing(100,30);
    cst{getStructIdx('L Eye'),6}{1}        = DoseObjectives.matRad_SquaredOverdosing(100,30);
    cst{getStructIdx('patient'),6}{1}      = DoseObjectives.matRad_SquaredOverdosing(10,40);

    pln.radiationMode   = 'protons';
    pln.machine         = 'Generic';
    pln.bioModel        = 'none';
    pln.multScen        = 'nomScen';
    pln.numOfFractions  = 30;

    pln.propStf.gantryAngles = [90 270];
    pln.propStf.couchAngles  = zeros(size(pln.propStf.gantryAngles));
    pln.propStf.bixelWidth   = 3;
    pln.propStf.isoCenter    = matRad_getIsoCenter(cst,ct,0);

    pln.propDoseCalc.engine = 'HongPB';
    pln.propDoseCalc.doseGrid.resolution.x = 3;
    pln.propDoseCalc.doseGrid.resolution.y = 3;
    pln.propDoseCalc.doseGrid.resolution.z = 3;

    pln.propOpt.quantityOpt = 'physicalDose';
    pln.propOpt.optimizer   = 'IPOPT';
    pln.propOpt.runDAO      = false;

    stf       = matRad_generateStf(ct,cst,pln);
    dij       = matRad_calcDoseInfluence(ct,cst,stf,pln);
    resultGUI = matRad_fluenceOptimization(dij,cst,pln);
    clear dij;

    doseHong = resultGUI.physicalDose;
    w        = resultGUI.w;

    save(planFile,'ct','cst','pln','stf','w','doseHong','numHistoriesDirect');
else
    load(planFile,'ct','cst','pln','stf','w','doseHong');
end

%% TOPAS settings
plnTopas = pln;
% the pln still selects the Hong engine used for the optimization
plnTopas.propDoseCalc = rmfield(plnTopas.propDoseCalc,'engine');

% same material converter as the other TOPAS recalculations
plnTopas.propDoseCalc.materialConverter.loadConverterFromFile    = true;
plnTopas.propDoseCalc.infilenames.matConv_Schneider_loadFromFile = 'TOPAS_SchneiderConverter_MCsquareDensity.txt.in';

%% TOPAS: write input files
if writeFiles
    engineTopasWrite = DoseEngines.matRad_TopasMCEngine(plnTopas);
    engineTopasWrite.numHistoriesDirect    = numHistoriesDirect;
    engineTopasWrite.externalCalculation   = 'write';
    engineTopasWrite.scorer.reportQuantity = {'Sum'};

    resultTopasWrite = engineTopasWrite.calcDoseForward(ct,cst,stf,w);
    movefile(fullfile(resultTopasWrite.meta.TOPASworkingDir,'*'),topasDir);
end

%% TOPAS: run
if runTopas
    topasPlanFiles = dir(fullfile(topasDir,'matRad_plan_field*_run*.txt'));
    if isempty(topasPlanFiles)
        matRad_cfg.dispError('No TOPAS plan files (matRad_plan_field*_run*.txt) found in %s',topasDir);
    end

    startDir = pwd;
    cd(topasDir);
    for i = 1:numel(topasPlanFiles)
        [~,fname] = fileparts(topasPlanFiles(i).name);
        topasCall = sprintf('%s %s.txt > %s.out > %s.log',topasExecCommand,fname,fname,fname);
        status = system(topasCall,'-echo');
        if status ~= 0
            cd(startDir);
            matRad_cfg.dispError('TOPAS exited with error code %d for %s',status,fname);
        end
    end
    cd(startDir);
end

%% TOPAS: read results
if readTopas
    engineTopas = DoseEngines.matRad_TopasMCEngine(plnTopas);
    engineTopas.externalCalculation = topasDir;
    resultTopas = engineTopas.calcDoseForward(ct,cst,stf,w);
    doseTopas   = resultTopas.physicalDose;

    save(fullfile(mainDir,'protons_TG166brain_TOPAS.mat'),'ct','cst','pln','stf','w','doseHong','doseTopas');

    %% Compare HongPB and TOPAS: dose slices and difference
    isoIdx     = matRad_world2cubeIndex(pln.propStf.isoCenter,ct); % [yIdx xIdx zIdx]
    doseWindow = [0 max([doseHong(:); doseTopas(:)])];
    doseDiff   = doseTopas - doseHong;
    diffWindow = [-1 1]*max(abs(doseDiff(:)));

    figure('Color','w','WindowState','maximized');
    subplot(1,3,1);
    matRad_plotSlice(ct,'axesHandle',gca,'cst',cst,'dose',doseHong,'plane',3,'slice',isoIdx(3), ...
        'doseWindow',doseWindow,'title','HongPB');
    subplot(1,3,2);
    matRad_plotSlice(ct,'axesHandle',gca,'cst',cst,'dose',doseTopas,'plane',3,'slice',isoIdx(3), ...
        'doseWindow',doseWindow,'title','TOPAS');
    subplot(1,3,3);
    matRad_plotSlice(ct,'axesHandle',gca,'cst',cst,'dose',doseDiff,'plane',3,'slice',isoIdx(3), ...
        'doseWindow',diffWindow,'title','TOPAS - HongPB');

    %% Compare HongPB and TOPAS: DVHs
    dvhDoseGrid = linspace(0,1.05*doseWindow(2),1000);
    dvhHong  = matRad_calcDVH(cst,doseHong,'cum',dvhDoseGrid);
    dvhTopas = matRad_calcDVH(cst,doseTopas,'cum',dvhDoseGrid);

    colors = lines(size(cst,1));
    figure('Color','w');
    hold on;
    for i = 1:size(cst,1)
        plot(dvhDoseGrid,dvhHong(i).volumePoints,'-','Color',colors(i,:),'LineWidth',1.5, ...
            'DisplayName',cst{i,2});
        plot(dvhDoseGrid,dvhTopas(i).volumePoints,'--','Color',colors(i,:),'LineWidth',1.5, ...
            'HandleVisibility','off');
    end
    grid on; box on;
    ylim([0 100]);
    xlabel('physical dose [Gy]');
    ylabel('volume [%]');
    title('DVH: HongPB (solid) vs. TOPAS (dashed)');
    legend('show','Location','northeastoutside');
end
