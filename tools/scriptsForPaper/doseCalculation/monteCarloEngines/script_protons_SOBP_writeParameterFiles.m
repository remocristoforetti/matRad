%% Script: Write TOPAS/MCsquare/FRED input files for the optimized SOBP plan
%% set matRad runtime configuration
clear;
matRad_rc;
matRad_cfg = MatRad_Config.instance();

%% Load the SOBP plan optimized in analyticalEngines/script_protons_SOBP_analytical.m
load(fullfile(matRad_cfg.primaryUserFolder,'dataPaper','doseCalculation','analyticalEngines','protons_SOBP_analytical_2.mat'), ...
    'ct','cst','pln','stf','w');
mainDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','doseCalculation','MCengines','protonsSOBP_2');
if ~isfolder(mainDir)
    mkdir(mainDir);
end

% the pln still selects the Hong engine used for the optimization
if isfield(pln.propDoseCalc,'engine')
    pln.propDoseCalc = rmfield(pln.propDoseCalc,'engine');
end

numHistoriesDirect = 1e6;

%% Save the plan information for the run/read scripts
planFile = fullfile(mainDir,'planInfo.mat');
save(planFile,'ct','cst','pln','stf','w','numHistoriesDirect');

%% TOPAS: write input files
topasWorkingDir = fullfile(mainDir,'TOPAS');
if ~isfolder(topasWorkingDir)
    mkdir(topasWorkingDir);
end

pln.propDoseCalc.materialConverter.loadConverterFromFile     = true;
pln.propDoseCalc.infilenames.matConv_Schneider_loadFromFile  = 'TOPAS_SchneiderConverter_MCsquareDensity.txt.in';

engineTopasWrite = DoseEngines.matRad_TopasMCEngine(pln);
engineTopasWrite.numHistoriesDirect = numHistoriesDirect;
engineTopasWrite.externalCalculation = 'write';
engineTopasWrite.scorer.reportQuantity = {'Sum'};

resultTopasWrite = engineTopasWrite.calcDoseForward(ct,cst,stf,w);
movefile(fullfile(resultTopasWrite.meta.TOPASworkingDir,'*'),topasWorkingDir);

%% MCsquare: write input files
mcSquareWorkingDir = fullfile(mainDir,'MCsquare');
if ~isfolder(mcSquareWorkingDir)
    mkdir(mcSquareWorkingDir);
end

engineMCsquareWrite = DoseEngines.matRad_ParticleMCsquareEngine(pln);
engineMCsquareWrite.workingDir = mcSquareWorkingDir;
engineMCsquareWrite.numHistoriesDirect = numHistoriesDirect;
engineMCsquareWrite.externalCalculation = 'write';
engineMCsquareWrite.calcDoseForward(ct,cst,stf,w);

%% FRED: write input files
fredWorkingDir = fullfile(mainDir,'FRED');
if ~isfolder(fredWorkingDir)
    mkdir(fredWorkingDir);
end

engineFredWrite = DoseEngines.matRad_ParticleFREDEngine(pln);
engineFredWrite.HUtable = 'matRad_MCsquareDensity_FredMaterialConverter';
engineFredWrite.workingDir = fredWorkingDir;
engineFredWrite.numHistoriesDirect = numHistoriesDirect;
engineFredWrite.externalCalculation = 'write';

engineFredWrite.primaryMass = 1;
engineFredWrite.calcDoseForward(ct,cst,stf,w);
