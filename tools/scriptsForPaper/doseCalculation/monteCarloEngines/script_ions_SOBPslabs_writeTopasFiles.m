%% Script: Write TOPAS input files for the optimized helium/carbon/oxygen SOBP plans (entrance slabs)
%% set matRad runtime configuration
clear;
matRad_rc;
matRad_cfg = MatRad_Config.instance();

ions = {'helium','carbon','oxygen'};
numHistoriesDirect = 1e6;

for i = 1:numel(ions)
    %% Load the plan optimized in analyticalEngines/script_ions_SOBPslabs_analytical.m
    load(fullfile(matRad_cfg.primaryUserFolder,'dataPaper','doseCalculation','analyticalEngines',[ions{i} '_SOBPslabs_analytical_2.mat']), ...
        'ct','cst','pln','stf','w');
    mainDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','doseCalculation','MCengines',[ions{i} 'SOBPslabs_2']);
    if ~isfolder(mainDir)
        mkdir(mainDir);
    end

    % the pln still selects the Hong engine used for the optimization
    if isfield(pln.propDoseCalc,'engine')
        pln.propDoseCalc = rmfield(pln.propDoseCalc,'engine');
    end

    %% Save the plan information for the run/read scripts
    save(fullfile(mainDir,'planInfo.mat'),'ct','cst','pln','stf','w','numHistoriesDirect');

    %% TOPAS: write input files
    topasWorkingDir = fullfile(mainDir,'TOPAS');
    if ~isfolder(topasWorkingDir)
        mkdir(topasWorkingDir);
    end

    engineTopasWrite = DoseEngines.matRad_TopasMCEngine(pln);
    engineTopasWrite.numHistoriesDirect = numHistoriesDirect;
    engineTopasWrite.externalCalculation = 'write';
    engineTopasWrite.scorer.reportQuantity = {'Sum'};

    resultTopasWrite = engineTopasWrite.calcDoseForward(ct,cst,stf,w);
    movefile(fullfile(resultTopasWrite.meta.TOPASworkingDir,'*'),topasWorkingDir);
end
