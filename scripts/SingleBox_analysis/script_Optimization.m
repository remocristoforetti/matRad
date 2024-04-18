%% load parameters
matRad_cfg = MatRad_Config.instance();
load('BOXPHANTOM_OAR_dis.mat');
saveDir = fullfile(matRad_cfg.matRadRoot, 'BOXPHANTOM_OAR_dist_analysis', '2mm');


load(fullfile(saveDir, 'probQuantities.mat'),'exp', 'omega','dij', 'pln', 'cst');
load(fullfile(saveDir, 'stf.mat'));

dij.physicalDoseOmega = omega;
dij.physicalDoseExp   = exp;
dij.physicalDose = {[]};
%% setup cost functions
script_OptimizationParameters_PROB_box;


cst{1,6}{1}.robustness = 'PROB';

cst{2,6}{1}.robustness = 'PROB';
cst{3,6}{1}.robustness = 'PROB';

Case = 7;
switch Case
    case 1
        cst{1,6} = {};
        cst{3,6} = {};

       
        cst{2,6}{2} = OmegaObjectives.matRad_TotalVariance();
        cst{2,6}{2}.penalty = cst{2,6}{1}.penalty;

    case 2
        
        cst{1,6}{2} = OmegaObjectives.matRad_TotalVariance();
        cst{1,6}{2}.penalty = cst{1,6}{1}.penalty;  
        
        cst{2,6}{2} = OmegaObjectives.matRad_TotalVariance();
        cst{2,6}{2}.penalty = cst{2,6}{1}.penalty;

        cst{3,6}{2} = OmegaObjectives.matRad_TotalVariance();
        cst{3,6}{2}.penalty = cst{3,6}{1}.penalty;

    case 5

        cst{2,6}{2} = OmegaConstraints.matRad_MaxVariance(0.0012);
        cst{2,6}{2}.robustness = 'PROB';

    case 7
        cst{2,6}{2} = OmegaConstraints.matRad_MaxVariance(0.0010);
        cst{2,6}{2}.robustness = 'PROB';
end

%% Optimization
matRad_cfg = MatRad_Config.instance();
matRad_cfg.propOpt.defaultMaxIter = 3000;
pln.propOpt.visualizeSingleObjectives = false;

resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% save
w = resultGUI.w;
prob_info = resultGUI.info;

costFunctions = resultGUI.costFunctions;
save(fullfile(saveDir,['planResults_', num2str(Case), '.mat']), 'cst', 'w', 'prob_info','costFunctions', '-v7.3');