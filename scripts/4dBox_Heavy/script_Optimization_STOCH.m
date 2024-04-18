%% load parameters
matRad_cfg = MatRad_Config.instance();
load('4D_BOXPHANTOM_OAR_dist_heavy.mat');

saveDir = fullfile(matRad_cfg.matRadRoot, '4D_BOXPHANTOM_OAR_dist_heavy_analysis', '2mm');


load(fullfile(saveDir, 'probQuantities.mat'), 'pln');
% load(fullfile(saveDir, 'stf.mat'));
% 
% dij.physicalDoseOmega = omega;
% dij.physicalDoseExp   = exp;
% dij.physicalDose = {[]};
dij = matRad_loadDijScenarios(ct,saveDir,'all');
%% setup cost functions
script_OptimizationParameters_PROB_4Dbox;


cst{1,6}{1}.robustness = 'STOCH';

cst{2,6}{1}.robustness = 'STOCH';
cst{3,6}{1}.robustness = 'STOCH';

Case = 4;
switch Case
    case 3
        cst{1,6} = {};
        cst{3,6} = {};

    case 4
        
end

%% Optimization
matRad_cfg = MatRad_Config.instance();
matRad_cfg.propOpt.defaultMaxIter = 1000;
pln.propOpt.scen4D = 'all';
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% save
w = resultGUI.w;
prob_info = resultGUI.info;

costFunctions = resultGUI.costFunctions;
save(fullfile(saveDir,['planResults_', num2str(Case), '.mat']), 'cst', 'w', 'prob_info','costFunctions', '-v7.3');