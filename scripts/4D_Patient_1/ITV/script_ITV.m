%% load ct
matRad_cfg = MatRad_Config.instance();
load('C:\r408i_data\r408i_data\CTDatasetMotion\102_HM10395_333_newCst.mat');

matRad_cfg = MatRad_Config.instance();

saveDir = fullfile(matRad_cfg.matRadRoot, 'Patient_1_analysis', 'ITV');
%% Add ITV
[cst,ctITV] = matRad_calcITV(cst,ct);


%% ADD margin
itvMargin.x = 4; % mm
itvMargin.y = 4; % mm
itvMargin.z = 4; % mm

itvMask = zeros(ct.cubeDim);
itvMask(cst{14,4}{1}) = 1;
expItvMask = matRad_addMargin(itvMask,cst,ct.resolution,itvMargin);


cst{15,1} = 14;
cst{15,2} = 'ITVmargin';
cst{15,3} = 'TARGET';
cst{15,4} = repmat({find(expItvMask)}, 1, ct.numOfCtScen);
cst{15,5} = struct('Priority', 1, ...
                    'alphaX', 0.1, ...
                    'betaX', 0.05, ...
                    'Visible', 1, ...
                    'visibleColor', [1 1 1]);

cst{15,6} = cst{8,6};
cst{8,6} = [];
%% Plan setup
script_planParameters_Patient;
pln.multScen = matRad_multScen(ct, 'nomScen');

pln.propDoseCalc.doseGrid.resolution.x = 2;
pln.propDoseCalc.doseGrid.resolution.y = 2;
pln.propDoseCalc.doseGrid.resolution.z = 2;

%% stf
stf = matRad_generateStf(ctITV,cst,pln);
%% Dose calculation
dij = matRad_calcParticleDoseMultipleScenarios(ct, stf,pln,cst,0,saveDir,0);

%% Load scenarios
dij = matRad_loadDijScenarios(ct,saveDir,'all');

%% cst
script_OptimizationParameters_PROB_Patient;

cst{3,6}{1}.robustness = 'none';
cst{4,6}{1}.robustness = 'none';

Case = 200;

switch Case

    case 20

        cst{15,6} = cst{8,6};

        cst{8,6} = [];

        %cst{13,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(500,18));
        %cst{13,6}{1}.robustness = 'none';
    case 21

        cst{15,6}{1} = struct(DoseObjectives.matRad_MeanDose(100,15));

        cst{8,6} = [];

    case 200
        cst{3,6} = {};
        cst{4,6} = {};
        cst{8,6} = {};
        cst{13,6} = {};

        %cst{3,6}{1} = struct(DoseObjectives.matRad_MeanDose(100,20, 'Quadratic'));

        cst{3,6}{2} = struct(DoseConstraints.matRad_MinMaxMeanDose(0,20));
        cst{3,6}{2}.robustness = 'none';
        cst{3,6}{1} = struct(DoseObjectives.matRad_MaxDVH(100,30,46));
        cst{3,6}{3} = struct(DoseObjectives.matRad_MeanDose(10,0));
        cst{3,6}{3}.robustness = 'none';
        %cst{3,6}{1}.robustness = 'PROB';
        cst{3,6}{1}.robustness = 'none';
   


        %cst{4,6}{1} = struct(DoseObjectives.matRad_MeanDose(100,15, 'Quadratic'));
        cst{4,6}{2} = struct(DoseConstraints.matRad_MinMaxMeanDose(0,15));
        cst{4,6}{2}.robustness = 'noen';
        cst{4,6}{1} = struct(DoseObjectives.matRad_MaxDVH(100,20,35));
        cst{4,6}{3} = struct(DoseObjectives.matRad_MeanDose(10,0));
        cst{4,6}{3}.robustness = 'none';

        %cst{4,6}{1}.robustness = 'PROB';
        cst{4,6}{1}.robustness = 'none';

        cst{15,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(100,60));
        cst{15,6}{1}.robustness = 'none';
end

%% Opt

pln.propOpt.scen4D = 1;

resultGUI = matRad_fluenceOptimization(dij,cst,pln);


%% Save
w = resultGUI.w;

prob_info = resultGUI.info;


costFunctions = resultGUI.costFunctions;

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

save(fullfile(saveDir,['ITV4mm_planResults_', num2str(Case), '.mat']), 'cst', 'w', 'prob_info','costFunctions','pln', 'stf', '-v7.3');