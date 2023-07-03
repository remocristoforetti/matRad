%% Define file with info
[fileName, folder] = uigetfile

%% Instantiate simulation class
doseSimulation = matRad_baseDataGeneration_dose;

doseSimulation.retriveMainClass();