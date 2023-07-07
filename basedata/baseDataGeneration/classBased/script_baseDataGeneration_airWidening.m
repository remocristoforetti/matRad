%% Instantiate the class and load parameters

matRad_cfg = MatRad_Config.instance();

[fileName,folder] = uigetfile([matRad_cfg.matRadRoot, filesep, 'baseData',filesep, 'baseDataGeneration','.mat'], 'Select Parameter file for main baseData');

%% Instantiate the class
baseData_airWidening = matRad_baseDataGeneration_airWidening();

baseData_airWidening.retriveMainClass([folder, fileName]);
%% Define scorers parameters
baseData_airWidening.scorerParams.scorers = {'PhaseSpace'};%{'DoseToMedium', 'Fluence','EdBinned', 'ProtonLET'};
baseData_airWidening.scorerParams.nScorers = size(baseData_airWidening.scorerParams.scorers,2);
baseData_airWidening.scorerParams.ions = {'protons'};
baseData_airWidening.scorerParams.ionsZ = 1;
% airWideningSimulation.scorerParams.Ebinning.EMin = 0;
% airWideningSimulation.scorerParams.Ebinning.EMax = 500;
% airWideningSimulation.scorerParams.Ebinning.nEBins = 1000;
%% Define phantoms
baseData_airWidening.phantoms.nPhantoms = 9;
baseData_airWidening.phantoms.depths    = [0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000];
baseData_airWidening.phantoms.HL        = 0.05; %mm
baseData_airWidening.phantoms.rMax      = 70; %mm

baseData_airWidening.phantoms.sourcePosition             = -baseData_airWidening.MCparams.BAMtoISO;
baseData_airWidening.MCparams.runDirectory = [baseData_airWidening.workingDir, filesep, 'SimulationAirWidening'];
baseData_airWidening.MCparams.nPrimaries = 250000;

baseData_airWidening.saveParameters();

%% Set fileName saved
fileName = [baseData_airWidening.workingDir, filesep, 'baseDataParameters', filesep,'AirWideningSimulation06-Jul-2023proton.mat'];
%% Instantiate the simulation subclass
% this should then be a subclass method that calls all the functions.
% Every subclass can in principle have different files to write
airWideningSimulation = matRad_airWidening_simulation();
airWideningSimulation.retriveMainClass(fileName);


airWideningSimulation.generateTreeDirectory();

airWideningSimulation.writeRunFiles();
airWideningSimulation.writeSimulationParameters();
airWideningSimulation.writeBasicFile();

airWideningSimulation.writeScorers();
%% save parameters
airWideningSimulation.saveParameters();

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Execute the analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Instantiate the class
airWideningAnalysis = matRad_airWidening_analysis();

airWideningAnalysis.retriveMainClass(fileName);

airWideningAnalysis.performAnalysis();

%% save the class object and output

airWideningAnalysis.saveOutput();