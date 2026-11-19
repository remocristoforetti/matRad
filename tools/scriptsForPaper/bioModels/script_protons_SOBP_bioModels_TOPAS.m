%% Script: TOPAS recalculation of the McNamara-optimized RBExDose SOBP
% Takes the plan optimized in script_protons_SOBP_bioModels_analytical.m
% and recalculates it with TOPAS using its RBE scorers (for MCN/WED this
% also scores LETd). Three steps, each toggled below: write the TOPAS
% input files, run TOPAS, read the results back into matRad.
% Note: the TOPAS engine does not use pln.bioModel, the model is selected
% with engine.scorer.RBE_model ('MCN', 'WED', 'LEM', 'libamtrack').

clear;
matRad_rc;
matRad_cfg = MatRad_Config.instance();

bioModel = 'MCN';

writeFiles = true;
runTopas   = true;
readTopas  = true;

numHistoriesDirect = 1e6;

topasExecCommand = 'wsl export TOPAS_G4_DATA_DIR=~/Applications/GEANT4/G4DATA; ~/Applications/TOPAS/OpenTOPAS-install/bin/topas';

%% Load the optimized plan
anaFile = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','bioModels','protons_SOBP_bioModels_analytical_2.mat');
load(anaFile,'ct','cst','pln','stf','w','physicalDose','LET','RBExDose');

pln.bioModel = bioModel;
% the pln still selects the Hong engine used for the optimization
if isfield(pln.propDoseCalc,'engine')
    pln.propDoseCalc = rmfield(pln.propDoseCalc,'engine');
end

% same material converter as the physical-dose TOPAS SOBP runs
pln.propDoseCalc.materialConverter.loadConverterFromFile    = true;
pln.propDoseCalc.infilenames.matConv_Schneider_loadFromFile = 'TOPAS_SchneiderConverter_MCsquareDensity.txt.in';

topasDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','bioModels','TOPAS',bioModel);
if ~isfolder(topasDir)
    mkdir(topasDir);
end

%% TOPAS: write input files
if writeFiles
    engineTopasWrite = DoseEngines.matRad_TopasMCEngine(pln);
    engineTopasWrite.numHistoriesDirect    = numHistoriesDirect;
    engineTopasWrite.externalCalculation   = 'write';
    engineTopasWrite.scorer.reportQuantity = {'Sum'};
    engineTopasWrite.scorer.RBE            = true;
    engineTopasWrite.scorer.RBE_model      = {bioModel};

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
    engineTopas = DoseEngines.matRad_TopasMCEngine(pln);
    engineTopas.externalCalculation = topasDir;
    resultTopas = engineTopas.calcDoseForward(ct,cst,stf,w);

    % RBE results are stored per TOPAS RBE model, e.g. RBExDose_MCN
    physicalDoseTopas = resultTopas.physicalDose;
    RBExDoseTopas     = resultTopas.(['RBExDose']);
    LETTopas          = resultTopas.LET .* resultTopas.physicalDose;
    LETDose           = LET .* physicalDose;

    save(fullfile(matRad_cfg.primaryUserFolder,'dataPaper','bioModels',sprintf('protons_SOBP_bioModels_TOPAS_%s.mat',bioModel)), ...
        'ct','cst','pln','stf','w','bioModel','physicalDoseTopas','RBExDoseTopas','LETTopas');

    %% Compare with the analytical (HongPB) calculation, same model
    isoIdx = [25,1,25]; %matRad_world2cubeIndex(pln.propStf.isoCenter,ct); % [yIdx xIdx zIdx]
    yIdx = isoIdx(1) + [-5:5];% + [-1,0,1];
    zIdx = isoIdx(3) + [-5:5];% + [-1,0,1];
    getProfile = @(cube) squeeze(mean(cube(yIdx,:,zIdx), [1,3]));

    figure('Color','w');
    yyaxis left;
    hold on;
    plot(ct.x, getProfile(physicalDose),            '-',  'LineWidth', 2, 'DisplayName', 'physical dose HongPB');
    plot(ct.x, getProfile(physicalDoseTopas),       '--', 'LineWidth', 2, 'DisplayName', 'physical dose TOPAS');
    plot(ct.x, getProfile(RBExDose.(bioModel)),     '-',  'LineWidth', 2, 'DisplayName', ['RBExDose HongPB ' bioModel]);
    plot(ct.x, getProfile(RBExDoseTopas),           '--', 'LineWidth', 2, 'DisplayName', ['RBExDose TOPAS ' bioModel]);
    ylabel('dose [Gy] / RBExDose [Gy(RBE)]');
    yyaxis right;
    plot(ct.x, getProfile(LETDose),  ':',  'LineWidth', 1.5, 'DisplayName', 'LET_d x D HongPB');
    plot(ct.x, getProfile(LETTopas), '-.', 'LineWidth', 1.5, 'DisplayName', 'LET_d x D TOPAS');
    ylabel('LET_d x D [keV/\mum Gy]');
    grid on;
    xlabel('x [mm]');
    legend('show','Location','northwest');
end
%% Tickz
