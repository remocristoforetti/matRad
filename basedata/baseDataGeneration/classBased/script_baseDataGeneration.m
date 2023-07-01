%% Instantiate the initial informaton required for the simulations
%load the machine as a starting point
load('protons_Generic.mat');
protons_machine = machine;
clear machine;

%% Instantiate the baseDataGenaration class for airWidening simulation
%airWideningSimulation = matRad_airWideningBaseData();
mainBaseDataGenerationClass = matRad_baseDataGeneration();
%% Specify energies to be simulated
simulateEnergies = [protons_machine.data(1:10:114).energy];
[~,eIdx] = intersect([protons_machine.data.energy], simulateEnergies);

%Load them into the class
mainBaseDataGenerationClass.simulateEnergies = simulateEnergies';

%% Define initFocus

initFocus = [];
%The starting point for beam widening is 
for k=1:size(simulateEnergies,2)
    initFocus.initSigma(k,1) = [protons_machine.data(eIdx(k)).initFocus.sigma(1)];
    initFocus.initThetaSigma(k,1) = 0;
    initFocus.correlation(k,1) = 0;
end

mainBaseDataGenerationClass.energyParams.initFocus = initFocus;
%% Define MC parameters
mainBaseDataGenerationClass.MCparams.runDirectory = [mainBaseDataGenerationClass.workingDir,filesep,'SimulationProtonsAirWidening'];
mainBaseDataGenerationClass.MCparams.nRuns = 4;               
mainBaseDataGenerationClass.MCparams.previousRuns = 0;     %additionaRuns to add (for later)
mainBaseDataGenerationClass.MCparams.doubleSource = 0;        %boolean to choose double source to add, for air widening is signe source
mainBaseDataGenerationClass.MCparams.nPrimaries = 250000;
mainBaseDataGenerationClass.MCparams.templateDir = [mainBaseDataGenerationClass.workingDir,filesep,'templates'];
mainBaseDataGenerationClass.MCparams.numberOfThreads = 0;
mainBaseDataGenerationClass.MCparams.BAMtoISO = 1000; %in mm
mainBaseDataGenerationClass.MCparams.sourceParticle = 'proton';

%% Save main class
mainBaseDataGenerationClass.saveParameters();