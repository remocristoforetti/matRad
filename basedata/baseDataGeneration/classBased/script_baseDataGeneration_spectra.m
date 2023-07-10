%% Instantiate the class and load parameters
matRad_cfg = MatRad_Config.instance();
[fileName,folder] = uigetfile([matRad_cfg.matRadRoot, filesep, 'baseData',filesep, 'baseDataGeneration','.mat'], 'Select Parameter file for main baseData');
%% Instantiate the class
spectraGeneration = matRad_baseDataGeneration_dose();

spectraGeneration.retriveMainClass([folder,fileName]);
%% Setup Phantom parameters
spectraGeneration.scorerParams.scorers = {'Fluence', 'EnergyDeposit'};%, 'Fluence','EdBinned', 'ProtonLET'};
spectraGeneration.scorerParams.nScorers = size(spectraGeneration.scorerParams.scorers,2);
spectraGeneration.scorerParams.ions = {'protons'};
spectraGeneration.scorerParams.ionsZ = 1;


%%% Phantoms could then become a class per se, so that will have a clas for
%%% the dose phantom, which is diveded in three and dorectly place it here
spectraGeneration.phantoms.nPhantoms = 3;

peakPhantomHL = 20*ones(spectraGeneration.energyParams.nEnergies,1); %mm
proximalPhantomHL = round((spectraGeneration.energyParams.simulateRanges-peakPhantomHL)/2);
proximalPhantomHL(proximalPhantomHL<1) = 1;
distalPhantomHL  = 25*ones(spectraGeneration.energyParams.nEnergies,1); %mm

spectraGeneration.phantoms.depths    = [proximalPhantomHL, (2*proximalPhantomHL + peakPhantomHL), 2*(proximalPhantomHL + peakPhantomHL)+distalPhantomHL];
spectraGeneration.phantoms.HL        = [proximalPhantomHL, peakPhantomHL, distalPhantomHL]; %mm
spectraGeneration.phantoms.rMax      = [50, 50, 50]; %mm
spectraGeneration.phantoms.Zbins     = [proximalPhantomHL,peakPhantomHL*10,distalPhantomHL]; % resolution [2 mm, 0.2 mm, 2 mm]
spectraGeneration.phantoms.Rbins     = [1, 1, 1];

spectraGeneration.phantoms.sourcePosition = -0.1; %mm
spectraGeneration.MCparams.runDirectory = [spectraGeneration.workingDir, filesep, 'SimulationSpectra'];
spectraGeneration.parameterVariableName = 'spectra';
%InitFocus information should be loaded here directly. This has to have
%same source position and parameters as a dose simulation
%% save parameters
spectraGeneration.saveParameters();


%% Setup simulation class

spectraSimulation = matRad_dose_simulation();

fileName = [spectraSimulation.workingDir, filesep, 'baseDataParameters', filesep, 'spectra05-Jul-2023proton'];
spectraSimulation.retriveMainClass(fileName);

%% Generate simulation files
spectraSimulation.generateTreeDirectory();
spectraSimulation.writeSimulationFiles();

%% Analysis
spectraAnalysis = matRad_spectra_analysis();
fileName = [spectraAnalysis.workingDir, filesep, 'baseDataParameters', filesep, 'spectra07-Jul-2023proton.mat'];
spectraAnalysis.retriveMainClass(fileName);

spectraAnalysis.MCparams.runDirectory = [spectraAnalysis.workingDir, filesep, 'SimulationSpectra'];

spectraAnalysis.performAnalysis();