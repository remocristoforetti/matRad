%% Define file with info
matRad_cfg = MatRad_Config.instance();

[fileName,folder] = uigetfile([matRad_cfg.matRadRoot, filesep, 'baseData',filesep, 'baseDataGeneration','.mat'], 'Select Parameter file for main baseData');


%% Instantiate simulation class
doseGeneration = matRad_baseDataGeneration_dose;

doseGeneration.retriveMainClass([folder, fileName]);

%% Instantiate subclass

doseSimulation = matRad_dose_simulation();


doseSimulation.retriveMainClass([folder,fileName]);

% Define MCparams
doseSimulation.MCparams.runDirectory = [doseSimulation.workingDir,filesep,'SimulationDose']; %-> This could go in the class constructor

%Retrive results for initFocus
initFocus = load([doseSimulation.workingDir,filesep, 'output', filesep,'airWideningAnalysis_04-Jul-2023_proton']);
initFocus = initFocus.saveStr.initFocus;


doseSimulation.interpInitFocus(initFocus);
doseSimulation.MCparams.doubleSource = 1;
%% Define phantoms and scorers
doseSimulation.scorerParams.scorers = {'DoseToMedium'};%, 'Fluence','EdBinned', 'ProtonLET'};
doseSimulation.scorerParams.nScorers = size(doseSimulation.scorerParams.scorers,2);
doseSimulation.scorerParams.ions = {'protons'};
doseSimulation.scorerParams.ionsZ = 1;


%%% Phantoms could then become a class per se, so that will have a clas for
%%% the dose phantom, which is diveded in three and dorectly place it here
doseSimulation.phantoms.nPhantoms = 3;

peakPhantomHL = 20*ones(doseSimulation.energyParams.nEnergies,1); %mm
proximalPhantomHL = round((doseSimulation.energyParams.simulateRanges-peakPhantomHL)/2);
proximalPhantomHL(proximalPhantomHL<1) = 1;
distalPhantomHL  = 25*ones(doseSimulation.energyParams.nEnergies,1); %mm

doseSimulation.phantoms.depths    = [proximalPhantomHL, (2*proximalPhantomHL + peakPhantomHL), 2*(proximalPhantomHL + peakPhantomHL)+distalPhantomHL];
doseSimulation.phantoms.HL        = [proximalPhantomHL, peakPhantomHL, distalPhantomHL]; %mm
doseSimulation.phantoms.rMax      = [50, 50, 50]; %mm
doseSimulation.phantoms.Zbins     = [proximalPhantomHL,peakPhantomHL*10,distalPhantomHL]; % resolution [2 mm, 0.2 mm, 2 mm]
doseSimulation.phantoms.Rbins     = [100, 100, 100];

doseSimulation.phantoms.sourcePosition = -0.1; %mm
%% Save parameters
doseSimulation.saveParameters();

%% write parameter files
doseSimulation.generateTreeDirectory();

doseSimulation.writeSimulationFiles();

%doseSimulation.writeRunFiles();
%doseSimulation.writeSimulationParameters();
%doseSimulation.writeBasicFile();

%doseSimulation.writeScorers();


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Execute the analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
doseAnalysis = matRad_dose_analysis();
fileName = [doseAnalysis.workingDir, filesep, 'baseDataParameters',filesep, 'doseSimulation04-Jul-2023proton'];
doseAnalysis.retriveMainClass(fileName);


%doseAnalysis.performAnalysis();
%% save output

doseAnalysis.saveOutput();

%% generate machine file
doseFit = load(['C:\r408i_data\r408i_data\matRad_varRBErobOpt\basedata\baseDataGeneration\output\', 'doseAnalysis_04-Jul-2023_proton.mat']);
doseFit = doseFit.saveStr;
airWideningFit = load(['C:\r408i_data\r408i_data\matRad_varRBErobOpt\basedata\baseDataGeneration\output\', 'airWideningAnalysis_04-Jul-2023_proton']);

airWideningFit = airWideningFit.saveStr;

machine = doseAnalysis.generateMachineFile(doseFit.fitDoseOutput,airWideningFit.initFocus);

%% Save machine
%save('protons_testClassGenericProton.mat', 'machine');
%% Visulaize PDDs
protonMachine = load('C:\r408i_data\r408i_data\matRad_varRBErobOpt\basedata\protons_Generic.mat');
protonMachine = protonMachine.machine;
[~, eIdx] = intersect([protonMachine.data(:).energy], doseSimulation.simulateEnergies);


figure;
for k=1:length(eIdx)
    plot(protonMachine.data(eIdx(k)).depths, protonMachine.data(eIdx(k)).Z, '.-', 'color', 'k');
    hold on;
    plot(machine.data(k).depths, machine.data(k).Z, '.-', 'color', 'r');
    
end


%% Vis lateral

figure;
for k=11:11%length(eIdx)
    plot(protonMachine.data(eIdx(k)).depths, protonMachine.data(eIdx(k)).sigma, '.-', 'color', 'k');
    hold on;
    plot(machine.data(k).depths, machine.data(k).sigma, '.-', 'color', 'r');
end