%% load ct
matRad_cfg = MatRad_Config.instance();
load('C:\r408i_data\r408i_data\CTDatasetMotion\102_HM10395_333_newCst.mat');

matRad_cfg = MatRad_Config.instance();

%% load plan parameters
script_planParameters_Patient;

%% Scenario setup
matRad_cfg = MatRad_Config.instance();

%% Compute scenarios
pln.multScen = matRad_multScen(ct,'rndScen');
pln.multScen.nSamples = 3;

stf = matRad_generateStf(ct,cst,pln);


saveDir = fullfile(matRad_cfg.matRadRoot, 'Patient_1_analysis', '3mm_30_patient');

[dij, dijTemplate]  = matRad_calcParticleDoseMultipleScenarios(ct,stf,pln,cst,0,saveDir,0);

%% save stf

load(fullfile(saveDir, 'stf.mat'), 'stf');

%% check scenarios
scen1 = load(fullfile(saveDir, 'scenario_1.mat'));
scen2 = load(fullfile(saveDir, 'scenario_3.mat'));

w = ones(dij.totalNumOfBixels,1);

scen1 = scen1.dijScenario{1}*w;
scen2 = scen2.dijScenario{1}*w;

