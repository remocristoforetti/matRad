%% load ct
matRad_cfg = MatRad_Config.instance();
load('C:\r408i_data\r408i_data\CTDatasetMotion\102_HM10395_333_newCst.mat');
saveDir = fullfile(matRad_cfg.matRadRoot, 'Patient_1_analysis', '3mm_30_patient');


matRad_cfg = MatRad_Config.instance();
%% load plan parameters
script_planParameters_Patient;

%% load scenarios

scenMode = 'none';
[dij, newMultiScen, loadedScenarios] = matRad_loadDijScenarios(ct, saveDir, scenMode);

pln.multScen = newMultiScen;
nScens = pln.multScen.totNumScen;


%% compute Omega
script_OptimizationParameters_PROB_Patient;
cstOnGrid = matRad_resizeCstToGrid(cst,dij.ctGrid.x, dij.ctGrid.y, dij.ctGrid.z, dij.doseGrid.x, dij.doseGrid.y, dij.doseGrid.z);

[exp, omega, omegaCalculationTime] = matRad_accumulateProbabilisticQuantities(saveDir,ct,cstOnGrid,pln,'all');


save(fullfile(saveDir,'probQuantities_30scens.mat'), 'exp', 'omega', 'dij', 'pln', 'cst', 'omegaCalculationTime', '-v7.3');