%% load parameters
matRad_cfg = MatRad_Config.instance();
load('C:\r408i_data\r408i_data\CTDatasetMotion\102_HM10395_333_newCst.mat');
saveDir = fullfile(matRad_cfg.matRadRoot, 'Patient_1_analysis', '2mm_30_patient');


script_planParameters_Patient;

[dij,newMultScen, ~] = matRad_loadDijScenarios(ct,saveDir,'all');
%% setup cost functions
script_OptimizationParameters_PROB_Patient;

cst{3,6}{1}.robustness = 'STOCH';
cst{4,6}{1}.robustness = 'STOCH';
cst{8,6}{1}.robustness = 'STOCH';

Case = 4;

switch Case
    case 3
       
        cst{1,6} = {};
        cst{3,6} = {};
        cst{4,6} = {};

    case 4
end

%% Optimization
matRad_cfg = MatRad_Config.instance();
matRad_cfg.propOpt.defaultMaxIter = 1000;
pln.propOpt.scen4D = 'all';
pln.multScen = newMultScen;


resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% save
w = resultGUI.w;
prob_info = resultGUI.info;

costFunctions = resultGUI.costFunctions;
save(fullfile(saveDir,['planResults_', num2str(Case), '.mat']), 'cst', 'w', 'prob_info','costFunctions','pln', '-v7.3');
