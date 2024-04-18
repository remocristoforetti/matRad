%% load parameters
matRad_cfg = MatRad_Config.instance();
load('4D_BOXPHANTOM_OAR_dist_heavy.mat');

% cubeIdx = [1,3,5];
% 
% cubes = [];
% cubesHU = [];
% 
% for i=1:3
% 
%     cubes{i}   = ct.cube{cubeIdx(i)};
%     cubesHU{i} = ct.cubeHU{cubeIdx(i)};
% end
% 
% ct.cube = cubes;
% ct.cubeHU = cubesHU;
% 
% ct.numOfCtScen = 3;

saveDir = fullfile(matRad_cfg.matRadRoot, '4D_BOXPHANTOM_OAR_dist_heavy_analysis', '2mm');


load(fullfile(saveDir, 'probQuantities.mat'),'exp', 'omega','dij', 'pln', 'cst');
load(fullfile(saveDir, 'stf.mat'));

dij.physicalDoseOmega = omega;
dij.physicalDoseExp   = exp;
dij.physicalDose = {[]};
%% setup cost functions
script_OptimizationParameters_PROB_4Dbox;


cst{1,6}{1}.robustness = 'PROB';

cst{2,6}{1}.robustness = 'PROB';
cst{3,6}{1}.robustness = 'PROB';

Case = 7;
switch Case
    case 1
        cst{1,6} = {};
        
        cst{3,6} = {};
        % cst{1,6}{2} = OmegaObjectives.matRad_TotalVariance();
        % cst{1,6}{2}.penalty = cst{1,6}{1}.penalty;  
        cst{2,6}{2} = OmegaObjectives.matRad_TotalVariance();
        cst{2,6}{2}.penalty = cst{2,6}{1}.penalty;

        % cst{3,6}{2} = OmegaObjectives.matRad_TotalVariance();
        % cst{3,6}{2}.penalty = cst{3,6}{1}.penalty;
    case 2
        
        cst{1,6}{2} = OmegaObjectives.matRad_TotalVariance();
        cst{1,6}{2}.penalty = cst{1,6}{1}.penalty;  
        
        cst{2,6}{2} = OmegaObjectives.matRad_TotalVariance();
        cst{2,6}{2}.penalty = cst{2,6}{1}.penalty;

        cst{3,6}{2} = OmegaObjectives.matRad_TotalVariance();
        cst{3,6}{2}.penalty = cst{3,6}{1}.penalty;

    case 5

        cst{1,6}{2} = OmegaObjectives.matRad_TotalVariance();
        cst{1,6}{2}.penalty = cst{1,6}{1}.penalty;  
        
        cst{2,6}{2} = OmegaObjectives.matRad_TotalVariance();
        cst{2,6}{2}.penalty = cst{2,6}{1}.penalty;

        cst{3,6}{2} = OmegaObjectives.matRad_TotalVariance();
        cst{3,6}{2}.penalty = 10*cst{3,6}{1}.penalty;

    case 7
        cst{1,6}{2} = OmegaObjectives.matRad_TotalVariance();
        cst{1,6}{2}.penalty = cst{1,6}{1}.penalty;  
        
        cst{2,6}{2} = OmegaObjectives.matRad_TotalVariance();
        cst{2,6}{2}.penalty = cst{2,6}{1}.penalty;

        cst{3,6}{2} = OmegaObjectives.matRad_TotalVariance();
        cst{3,6}{2}.penalty = 0.3*cst{3,6}{1}.penalty;
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
if any(size(dij.physicalDoseExp)>1)
    save(fullfile(saveDir,['planResults_phase_', num2str(Case), '.mat']), 'cst', 'w', 'prob_info','costFunctions', '-v7.3');
else
    save(fullfile(saveDir,['planResults_', num2str(Case), '.mat']), 'cst', 'w', 'prob_info','costFunctions', '-v7.3');
end