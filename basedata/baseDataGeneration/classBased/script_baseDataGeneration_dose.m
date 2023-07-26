%% Define file with info
matRad_cfg = MatRad_Config.instance();

[fileName,folder] = uigetfile([matRad_cfg.matRadRoot, filesep, 'baseData',filesep, 'baseDataGeneration','.mat'], 'Select Parameter file for main baseData');


%% Instantiate simulation class
doseGeneration = matRad_baseDataGeneration_dose;

doseGeneration.retriveMainClass([folder, fileName]);
%% Define phantoms and scorers
doseGeneration.scorers = {'DoseToMedium'};%, 'Fluence','EdBinned', 'ProtonLET'};
doseGeneration.scorerParams.nScorers = size(doseGeneration.scorers,2);
doseGeneration.scorerParams.ions = {'protons'};
doseGeneration.scorerParams.ionsZ = 1;

doseGeneration.scorerParams.filteredScorer = 0;

doseGeneration.scorerParams.energyBinned = 0;
%%% Phantoms could then become a class per se, so that will have a clas for
%%% the dose phantom, which is diveded in three and dorectly place it here
doseGeneration.phantoms.nPhantoms = 3;

peakPhantomHL = 20*ones(doseGeneration.energyParams.nEnergies,1); %mm
proximalPhantomHL = round((doseGeneration.energyParams.simulateRanges-peakPhantomHL)/2);
proximalPhantomHL(proximalPhantomHL<1) = 1;
distalPhantomHL  = 25*ones(doseGeneration.energyParams.nEnergies,1); %mm

doseGeneration.phantoms.material = 'MyWater';
doseGeneration.phantoms.depths    = [proximalPhantomHL, (2*proximalPhantomHL + peakPhantomHL), 2*(proximalPhantomHL + peakPhantomHL)+distalPhantomHL];
doseGeneration.phantoms.HL        = [proximalPhantomHL, peakPhantomHL, distalPhantomHL]; %mm
doseGeneration.phantoms.rMax      = [50, 50, 50]; %mm
doseGeneration.phantoms.Zbins     = [proximalPhantomHL,peakPhantomHL*10,distalPhantomHL]; % resolution [2 mm, 0.2 mm, 2 mm]
doseGeneration.phantoms.Rbins     = [100, 100, 100];

doseGeneration.phantoms.sourcePosition = -0.1; %mm

doseGeneration.parameterVariableName = 'doseGeneration';
doseGeneration.MCparams.runDirectory = [doseGeneration.workingDir, filesep, 'SimulationDose'];
%% Save parameters

doseGeneration.saveParameters();

%% Instantiate subclass
doseSimulation = matRad_dose_simulation();

fileName = [doseSimulation.workingDir, filesep, 'baseDataParameters', filesep, 'doseGeneration10-Jul-2023proton.mat'];

doseSimulation.retriveMainClass(fileName);
doseSimulation.MCparams.runDirectory = [matRad_cfg.matRadRoot, filesep,'baseData', filesep, 'baseDataGeneration', filesep, 'SimulationDose'];
% Define MCparams
doseSimulation.MCparams.runDirectory = [doseSimulation.workingDir,filesep,'SimulationDose']; %-> This could go in the class constructor

%Retrive results for initFocus
initFocus = load([doseSimulation.workingDir,filesep, 'output', filesep,'airWideningAnalysis_07-Jul-2023_proton.mat']);
initFocus = initFocus.saveStr.initFocus;


doseSimulation.interpInitFocus(initFocus);
doseSimulation.MCparams.doubleSource = 0;


%% write parameter files
doseSimulation.generateTreeDirectory();

doseSimulation.writeSimulationFiles();

%doseSimulation.writeRunFiles();
%doseSimulation.writeSimulationParameters();
%doseSimulation.writeBasicFile();

%doseSimulation.writeScorers();
doseSimulation.parameterVariableName = 'doseSimulation';
doseSimulation.saveParameters();
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Execute the analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

doseAnalysis = matRad_dose_analysis();
fileName = [doseAnalysis.workingDir, filesep, 'baseDataParameters', filesep, 'doseGeneration10-Jul-2023proton.mat'];
doseAnalysis.retriveMainClass(fileName);
doseAnalysis.MCparams.runDirectory = [matRad_cfg.matRadRoot, filesep,'baseData', filesep, 'baseDataGeneration', filesep, 'SimulationDoseG4_WATER_DoseToMedium'];


doseAnalysis.performAnalysis();
%% save output
%doseAnalysis.saveOutput();

%% generate machine file

%doseFit = load(['C:\r408i_data\r408i_data\matRad_varRBErobOpt\basedata\baseDataGeneration\output\', 'doseAnalysis_11-Jul-2023_proton.mat']);
%doseFit = doseFit.saveStr;
%airWideningFit = load(['C:\r408i_data\r408i_data\matRad_varRBErobOpt\basedata\baseDataGeneration\output\', 'airWideningAnalysis_06-Jul-2023_proton.mat']);

%airWideningFit = airWideningFit.saveStr;

%machine = doseAnalysis.generateMachineFile(doseFit.fitDoseOutput,airWideningFit.initFocus);

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

%% Vis sigmas
n = 7;
depths = machine.data(n).depths;
sigma = machine.data(n).sigma;
sigma1 = machine.data(n).sigma1;
sigma2 = machine.data(n).sigma2;
weight = machine.data(n).weight;



figure;
% plot(depths,sigma, '.-');
% grid on;
% grid minor;
% ylim([0,12]);


figure;
subplot(1,3,1);
plot(depths,sigma1, '.-');
grid on;
grid minor;
ylim([0,12]);

subplot(1,3,2);
plot(depths,sigma2, '.-');
grid on;
grid minor;
ylim([0,100]);

subplot(1,3,3);
plot(depths,weight, '.-');
grid on;
grid minor;

ylim([0,1]);