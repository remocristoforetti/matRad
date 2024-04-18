%% load ct
matRad_cfg = MatRad_Config.instance();
load('C:\r408i_data\r408i_data\CTDatasetMotion\102_HM10395_333_newCst.mat');

matRad_cfg = MatRad_Config.instance();

saveDir = fullfile(matRad_cfg.matRadRoot, 'Patient_1_analysis', 'ITV');
%% Plan setup
script_planParameters_Patient;

pln.propDoseCalc.doseGrid.resolution.x = 2;
pln.propDoseCalc.doseGrid.resolution.y = 2;
pln.propDoseCalc.doseGrid.resolution.z = 2;

load(fullfile(saveDir,['ITV4mm_planResults_', num2str(Case), '.mat']), 'cst', 'stf');

%% Dose calc

originalSaveDir = fullfile(matRad_cfg.matRadRoot, 'Patient_1_analysis', '2mm_30_patient');
originalCT = load('C:\r408i_data\r408i_data\CTDatasetMotion\102_HM10395_333_newCst.mat', 'ct');

originalCT = originalCT.ct;

[~,newMultiscen,~] = matRad_loadDijScenarios(originalCT, originalSaveDir, 'none');
pln.multScen = newMultiscen;

[dij, dijTemplate]  = matRad_calcParticleDoseMultipleScenarios(originalCT,stf,pln,cst,0,fullfile(saveDir, 'scenarios'),0);


%% See margin
[marginITV] = matRad_calcITVmarginForShiftScens(cstITV,ct, pln.multScen.isoShift.*3, 0.95);