%% Instantiate the initial informaton required for the simulations
%load the machine as a starting point
load('C:\r408i_data\r408i_data\JenniferScoringPipeline\Generic\protons_Generic.mat');
protons_machine = machine;
clear machine;

%% Instantiate the baseDataGenaration class for airWidening simulation
airWideningSimulation = matRad_airWideningBaseData();

%% Specify energies to be simulated
simulateEnergies = [protons_machine.data(1:10:114).energy];
[~,eIdx] = intersect([protons_machine.data.energy], simulateEnergies);

%Load them into the class
airWideningSimulation.simulateEnergies = simulateEnergies';

%% Define initFocus

initFocus = [];
%The starting point for beam widening is 
for k=1:size(simulateEnergies,2)
    initFocus.initSigma(k,1) = [protons_machine.data(eIdx(k)).initFocus.sigma(1)];
    initFocus.initThetaSigma(k,1) = 0;
    initFocus.correlation(k,1) = 0;
end

airWideningSimulation.energyParams.initFocus = initFocus;
%% Define MC parameters
airWideningSimulation.MCparams.runDirectory = [pwd,filesep,'SimulationProtonsAirWidening'];
airWideningSimulation.MCparams.nRuns = 4;               
airWideningSimulation.MCparams.previousRuns = 0;     %additionaRuns to add (for later)
airWideningSimulation.MCparams.doubleSource = 0;        %boolean to choose double source to add, for air widening is signe source
airWideningSimulation.MCparams.nPrimaries = 250000;
airWideningSimulation.MCparams.templateDir = [pwd,filesep,'templates'];
airWideningSimulation.MCparams.numberOfThreads = 0;
airWideningSimulation.MCparams.BAMtoISO = 1000; %in mm
airWideningSimulation.MCparams.sourceParticle = 'proton';
%% Define scorers parameters
airWideningSimulation.scorerParams.scorers = {'PhaseSpace'};%{'DoseToMedium', 'Fluence','EdBinned', 'ProtonLET'};
airWideningSimulation.scorerParams.nScorers = size(airWideningSimulation.scorerParams.scorers,2);
airWideningSimulation.scorerParams.ions = {'protons'};
airWideningSimulation.scorerParams.ionsZ = 1;
% airWideningSimulation.scorerParams.Ebinning.EMin = 0;
% airWideningSimulation.scorerParams.Ebinning.EMax = 500;
% airWideningSimulation.scorerParams.Ebinning.nEBins = 1000;
%% Define phantoms
airWideningSimulation.phantoms.nPhantoms = 5;
airWideningSimulation.phantoms.depths    = [0, 500, 1000, 1500, 2000];
airWideningSimulation.phantoms.HL        = 0.05; %mm
airWideningSimulation.phantoms.rMax      = 50; %mm
%% Start writing files
airWideningSimulation.generateTreeDirectory();

airWideningSimulation.writeRunFiles();
airWideningSimulation.writeSimulationParameters();
airWideningSimulation.writeBasicFile();

airWideningSimulation.writeScorers();

%% save
%airWideningSimulation.saveParameters();

%% Execute the analysis

airWideningAnalysis = matRad_baseDataAnalysisAirWidening(airWideningSimulation);

airWideningAnalysis.performAnalysis();

airWideningAnalysis.saveOutput();