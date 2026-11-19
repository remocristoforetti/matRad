%% Script: TOPAS recalculation of the McNamara-optimized TG166 brain plan
clear;
matRad_rc;
matRad_cfg = MatRad_Config.instance();

bioModel = 'MCN';

writeFiles = true;
runTopas   = false;
readTopas  = true;

numHistoriesDirect = 1e7;
numOfRuns          = 10; % histories are split over this many independent runs (one seed each) per field

topasExecCommand = 'wsl export TOPAS_G4_DATA_DIR=~/Applications/GEANT4/G4DATA; ~/Applications/TOPAS/OpenTOPAS-install/bin/topas';

%% Load the optimized plan
anaFile = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','bioModels','protons_TG166brain_MCN.mat');
load(anaFile,'ct','cst','pln','stf','w','physicalDose','LET','RBExDose');

pln.bioModel = bioModel;
% the pln still selects the Hong engine used for the optimization
if isfield(pln.propDoseCalc,'engine')
    pln.propDoseCalc = rmfield(pln.propDoseCalc,'engine');
end

% same material converter as the other TOPAS recalculations
pln.propDoseCalc.materialConverter.loadConverterFromFile    = true;
pln.propDoseCalc.infilenames.matConv_Schneider_loadFromFile = 'TOPAS_SchneiderConverter_MCsquareDensity.txt.in';

topasDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','bioModels','TOPAS',['TG166brain_' bioModel]);
if ~isfolder(topasDir)
    mkdir(topasDir);
end

%% TOPAS: write input files
if writeFiles
    engineTopasWrite = DoseEngines.matRad_TopasMCEngine(pln);
    engineTopasWrite.numHistoriesDirect    = numHistoriesDirect;
    engineTopasWrite.numOfRuns             = numOfRuns;
    engineTopasWrite.externalCalculation   = 'write';
    engineTopasWrite.scorer.reportQuantity = {'Sum'};
    engineTopasWrite.scorer.RBE            = true;
    engineTopasWrite.scorer.RBE_model      = {bioModel};
    engineTopasWrite.ignoreOutsideDensities = true;

    resultTopasWrite = engineTopasWrite.calcDoseForward(ct,cst,stf,w);
    movefile(fullfile(resultTopasWrite.meta.TOPASworkingDir,'*'),topasDir);

    % Cluster helpers: one line per TOPAS run (field x run) in topasRunList.txt,
    % and a SLURM array script running line SLURM_ARRAY_TASK_ID of that list.
    % Adapt the TOPAS environment lines in submitTopas.sh to the cluster.
    topasPlanFiles = dir(fullfile(topasDir,'matRad_plan_field*_run*.txt'));
    fid = fopen(fullfile(topasDir,'topasRunList.txt'),'w');
    fprintf(fid,'%s\n',topasPlanFiles.name);
    fclose(fid);

    % fid = fopen(fullfile(topasDir,'submitTopas.sh'),'w');
    % fprintf(fid,'#!/bin/bash\n');
    % fprintf(fid,'#SBATCH --job-name=topas_TG166_%s\n',bioModel);
    % fprintf(fid,'#SBATCH --array=1-%d\n',numel(topasPlanFiles));
    % fprintf(fid,'#SBATCH --cpus-per-task=1\n');
    % fprintf(fid,'#SBATCH --output=slurm_%%A_%%a.log\n\n');
    % fprintf(fid,'# export TOPAS_G4_DATA_DIR=/path/to/G4DATA\n');
    % fprintf(fid,'TOPAS=topas # /path/to/topas\n\n');
    % fprintf(fid,'cd "$(dirname "$0")"\n');
    % fprintf(fid,'planFile=$(sed -n "${SLURM_ARRAY_TASK_ID}p" topasRunList.txt)\n');
    % fprintf(fid,'$TOPAS "$planFile" > "${planFile%%.txt}.out"\n');
    % fclose(fid);

    % matRad_cfg.dispInfo('Wrote %d TOPAS runs (%d per field) to %s\n',numel(topasPlanFiles),numOfRuns,topasDir);
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

    engineTopas = DoseEngines.matRad_TopasMCEngine(pln);
    engineTopas.externalCalculation = topasDir;
    engineTopas.numOfRuns = numOfRuns;
    resultTopas = engineTopas.calcDoseForward(ct,cst,stf,w);

    physicalDoseTopas = resultTopas.physicalDose;
    RBExDoseTopas     = resultTopas.RBExDose;
    LETTopas          = resultTopas.LET;
    gamma             = matRad_gammaIndex(RBExDose, RBExDoseTopas, [ct.resolution.x, ct.resolution.y, ct.resolution.z], [2,2]);

    save(fullfile(matRad_cfg.primaryUserFolder,'dataPaper','bioModels',sprintf('protons_TG166brain_TOPAS_%s.mat',bioModel)), ...
        'ct','cst','pln','stf','w','bioModel','physicalDoseTopas','RBExDoseTopas','LETTopas');

    % Compare with the analytical (HongPB) calculation: slices
    isoIdx     = matRad_world2cubeIndex(pln.propStf.isoCenter,ct); % [yIdx xIdx zIdx]
    doseWindow = [0 max([RBExDose(:); RBExDoseTopas(:)])];
    letWindow  = [0 max(LET(:))];
    rbeDiff    = RBExDoseTopas - RBExDose;
    diffWindow = [-1 1]*max(abs(rbeDiff(:)));

    figure('Color','w','WindowState','maximized');
    subplot(2,3,1);
    matRad_plotSlice(ct,'axesHandle',gca,'cst',cst,'dose',RBExDose,'plane',3,'slice',isoIdx(3), ...
        'doseWindow',doseWindow,'title',['RBExDose HongPB ' bioModel]);
    subplot(2,3,2);
    matRad_plotSlice(ct,'axesHandle',gca,'cst',cst,'dose',RBExDoseTopas,'plane',3,'slice',isoIdx(3), ...
        'doseWindow',doseWindow,'title',['RBExDose TOPAS ' bioModel]);
    subplot(2,3,3);
    matRad_plotSlice(ct,'axesHandle',gca,'cst',cst,'dose',rbeDiff,'plane',3,'slice',isoIdx(3), ...
        'doseWindow',diffWindow,'title','RBExDose TOPAS - HongPB');
    subplot(2,3,4);
    matRad_plotSlice(ct,'axesHandle',gca,'cst',cst,'dose',LET,'plane',3,'slice',isoIdx(3), ...
        'doseWindow',letWindow,'title','LET_d HongPB [keV/\mum]');
    subplot(2,3,5);
    matRad_plotSlice(ct,'axesHandle',gca,'cst',cst,'dose',LETTopas,'plane',3,'slice',isoIdx(3), ...
        'doseWindow',letWindow,'title','LET_d TOPAS [keV/\mum]');
    subplot(2,3,6);
    matRad_plotSlice(ct,'axesHandle',gca,'cst',cst,'dose',gamma,'plane',3,'slice',isoIdx(3), ...
        'doseWindow',[0,2],'title','gamma', 'doseColorMap', gammaIndex);

    % Compare with the analytical (HongPB) calculation: DVHs
    dvhDoseGrid = linspace(0,1.05*doseWindow(2),1000);
    dvhHong  = matRad_calcDVH(cst,RBExDose,'cum',dvhDoseGrid);
    dvhTopas = matRad_calcDVH(cst,RBExDoseTopas,'cum',dvhDoseGrid);

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
    xlabel('RBExDose [Gy(RBE)]');
    ylabel('volume [%]');
    title(['RBExDose DVH ' bioModel ': HongPB (solid) vs. TOPAS (dashed)']);
    legend('show','Location','northeastoutside');
end
