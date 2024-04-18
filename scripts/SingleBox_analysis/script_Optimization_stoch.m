%% load parameters
matRad_cfg = MatRad_Config.instance();
load('BOXPHANTOM_OAR_dis.mat');
saveDir = fullfile(matRad_cfg.matRadRoot, 'BOXPHANTOM_OAR_dist_analysis', '2mm');


% load(fullfile(saveDir, 'probQuantities.mat'),'exp', 'omega','dij', 'pln', 'cst');
% load(fullfile(saveDir, 'stf.mat'));
% 
% dij.physicalDoseOmega = omega;
% dij.physicalDoseExp   = exp;
% dij.physicalDose = {[]};
dij = matRad_loadDijScenarios(ct,saveDir,'all');
%% setup cost functions
script_OptimizationParameters_PROB_box;

Case = 3;

switch Case
    case 3
       cst{1,6} = {};
       cst{3,6} = {};
       cst{2,6}{1}.robustness = 'STOCH';
       
    case 4

        cst{1,6}{1}.robustness = 'STOCH';
        cst{2,6}{1}.robustness = 'STOCH';
        cst{3,6}{1}.robustness = 'STOCH';
end


%% Optimization
matRad_cfg = MatRad_Config.instance();
matRad_cfg.propOpt.defaultMaxIter = 3000;

resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% save

w = resultGUI.w;
prob_info = resultGUI.info;

costFunctions = resultGUI.costFunctions;
save(fullfile(saveDir,['planResults_', num2str(Case), '.mat']), 'cst', 'w', 'prob_info','costFunctions', '-v7.3');