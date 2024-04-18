%% load parameters
matRad_cfg = MatRad_Config.instance();
load('C:\r408i_data\r408i_data\CTDatasetMotion\102_HM10395_333_newCst.mat');
%saveDir = fullfile(matRad_cfg.matRadRoot, 'ICCRdijs', '3mm_100');
% 
% 
% load(fullfile(saveDir, 'dij_Omega_100scens_3mm.mat'));
% load(fullfile(saveDir, 'stf.mat'));
% 
% dij = saveDij;
% dij.physicalDose = {[]};
% clear saveDij;
% script_planParameters_Patient;
% 
% [~,newMultScen, ~] = matRad_loadDijScenarios(ct,saveDir,'all', 'none');

saveDir = fullfile(matRad_cfg.matRadRoot, 'Patient_1_analysis', '2mm_30_patient');

[~,newMultScen, ~] = matRad_loadDijScenarios(ct,saveDir,'all', 'none');

load(fullfile(saveDir, 'probQuantities_30scens.mat'),'exp', 'omega','dij', 'pln', 'cst');
load(fullfile(saveDir, 'stf.mat'));

dij.physicalDoseOmega = omega;
dij.physicalDoseExp   = exp;
dij.physicalDose = {[]};
%% setup cost functions
script_OptimizationParameters_PROB_Patient;

%cst{1,6}{1}.robustness = 'PROB';
cst{3,6}{1}.robustness = 'PROB';
cst{4,6}{1}.robustness = 'PROB';
cst{8,6}{1}.robustness = 'PROB';
Case = 101;

switch Case
    case 1
       
        cst{1,6} = {};
        cst{3,6} = {};
        cst{4,6} = {};

        cst{8,6}{2} = OmegaObjectives.matRad_TotalVariance();
        cst{8,6}{2}.penalty = cst{8,6}{1}.penalty;

        % cst{3,6}{2} = OmegaObjectives.matRad_TotalVariance();
        % cst{3,6}{2}.penalty = cst{3,6}{1}.penalty;
    case 2
        
        %cst{1,6}{2} = OmegaObjectives.matRad_TotalVariance();
        %cst{1,6}{2}.penalty = cst{1,6}{1}.penalty;  
        
        cst{3,6}{2} = OmegaObjectives.matRad_TotalVariance();
        cst{3,6}{2}.penalty = cst{3,6}{1}.penalty;

        cst{4,6}{2} = OmegaObjectives.matRad_TotalVariance();
        cst{4,6}{2}.penalty = cst{4,6}{1}.penalty;

        cst{8,6}{2} = OmegaObjectives.matRad_TotalVariance();
        cst{8,6}{2}.penalty = cst{8,6}{1}.penalty;

    case 3
        cst{3,6}{2} = OmegaObjectives.matRad_TotalVariance();
        cst{3,6}{2}.penalty = cst{3,6}{1}.penalty;

        cst{4,6}{2} = OmegaObjectives.matRad_TotalVariance();
        cst{4,6}{2}.penalty = 10*cst{4,6}{1}.penalty;

        cst{8,6}{2} = OmegaObjectives.matRad_TotalVariance();
        cst{8,6}{2}.penalty = cst{8,6}{1}.penalty;
    case 4
        cst{3,6}{2} = OmegaConstraints.matRad_MaxVariance(0.01);
        cst{3,6}{2}.robustness = 'PROB';

        cst{4,6}{2} = OmegaConstraints.matRad_MaxVariance(0.01);
        cst{4,6}{2}.robustness = 'PROB';

    case 100

        cst{3,6} = {};
        cst{4,6} = {};
        cst{8,6} = {};
        cst{13,6} = {};

        %cst{3,6}{1} = struct(DoseObjectives.matRad_MeanDose(100,20, 'Quadratic'));

        cst{3,6}{2} = struct(DoseConstraints.matRad_MinMaxMeanDose(0,20));
        cst{3,6}{2}.robustness = 'PROB';
        cst{3,6}{1} = struct(DoseObjectives.matRad_MaxDVH(100,30,46));
        cst{3,6}{3} = struct(DoseObjectives.matRad_MeanDose(10,0));
        cst{3,6}{3}.robustness = 'PROB';
        %cst{3,6}{1}.robustness = 'PROB';
        cst{3,6}{1}.robustness = 'PROB';
   


        %cst{4,6}{1} = struct(DoseObjectives.matRad_MeanDose(100,15, 'Quadratic'));
        cst{4,6}{2} = struct(DoseConstraints.matRad_MinMaxMeanDose(0,15));
        cst{4,6}{2}.robustness = 'PROB';
        cst{4,6}{1} = struct(DoseObjectives.matRad_MaxDVH(100,20,35));
        cst{4,6}{3} = struct(DoseObjectives.matRad_MeanDose(10,0));
        cst{4,6}{3}.robustness = 'PROB';

        %cst{4,6}{1}.robustness = 'PROB';
        cst{4,6}{1}.robustness = 'PROB';

        cst{8,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(100,60));
        cst{8,6}{1}.robustness = 'PROB';
        cst{8,6}{2} = OmegaObjectives.matRad_TotalVariance();
        cst{8,6}{2}.robustness = 'PROB';
        cst{8,6}{2}.penalty = 100;

    case 101
                cst{3,6} = {};
        cst{4,6} = {};
        cst{8,6} = {};
        cst{13,6} = {};

        %cst{3,6}{1} = struct(DoseObjectives.matRad_MeanDose(100,20, 'Quadratic'));

        cst{3,6}{2} = struct(DoseConstraints.matRad_MinMaxMeanDose(0,20));
        cst{3,6}{2}.robustness = 'PROB';
        cst{3,6}{1} = struct(DoseObjectives.matRad_MaxDVH(100,30,46));
        cst{3,6}{3} = struct(DoseObjectives.matRad_MeanDose(10,0));
        cst{3,6}{3}.robustness = 'PROB';
        cst{3,6}{1}.robustness = 'PROB';
        cst{3,6}{4} = OmegaObjectives.matRad_TotalVariance();
        cst{3,6}{4}.penalty = 100;

   


        cst{4,6}{2} = struct(DoseConstraints.matRad_MinMaxMeanDose(0,15));
        cst{4,6}{2}.robustness = 'PROB';
        cst{4,6}{1} = struct(DoseObjectives.matRad_MaxDVH(100,20,35));
        cst{4,6}{3} = struct(DoseObjectives.matRad_MeanDose(10,0));
        cst{4,6}{3}.robustness = 'PROB';
        cst{4,6}{1}.robustness = 'PROB';

        %cst{4,6}{4} = OmegaObjectives.matRad_TotalVariance();
        %cst{4,6}{4}.penalty = 100;


        cst{8,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(100,60));
        cst{8,6}{1}.robustness = 'PROB';
        cst{8,6}{2} = OmegaObjectives.matRad_TotalVariance();
        cst{8,6}{2}.robustness = 'PROB';
        cst{8,6}{2}.penalty = 100;

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


