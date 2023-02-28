%% Setup Energies for LUT calculation
eMin = 0.1;
eMax = 300;
numOfEnergies = 100;
energies = logspace(log10(eMin),log10(eMax),numOfEnergies);

%% Setup the computation parameters:
track = Track_KC();
Z = [1 2 3 4 5 6];

domain = Domain();
domain.rNuc = 0.35;
domain.Beta_Tissue = 0.05;
domain.RN = 3.9;

%% Setup Homemade LUT calculator
InteGral = IntegralDose();
InteGral.Resolution = 0.01;
InteGral.stepWeightedIntegral = 0.01;
InteGral.corFact = 1;

%% Setup Survival class calculator
IntegralS = IntegralDose_Survival();

IntegralS.survivalSourcePath = '$HOME/Survival';

IntegralS.wDir = 'D:\matRad_gitHubRemo_MKM_TOPAS\MKM';
IntegralS.survivalParameterFileName = 'MKM_1';
IntegralS.calcProperties.projectName     = 'LUT_1';
IntegralS.calcProperties.output          = 'LQ_pars';
IntegralS.calcProperties.model           = 'MKM';
IntegralS.calcProperties.calculusType    = 'rapidMKM_Kase2008_corrected_beta';
%IntegralS.calcProperties.precision       = 0.15;
IntegralS.calcProperties.parallelismType = '0';
IntegralS.calcProperties.cellType        = 'HSG';

IntegralS.calcProperties.modelParam.alpha0   = 0.0;
IntegralS.calcProperties.modelParam.beta0    = 0.05;
IntegralS.calcProperties.modelParam.rNucleus = 3.1;
IntegralS.calcProperties.modelParam.rDomain  = 0.35;

IntegralS.calcProperties.ion = 'H';
IntegralS.calcProperties.trackMode = 'histogram';
IntegralS.calcProperties.energies  = energies;


idx = strfind(IntegralS.wDir, '\');
dirPath = IntegralS.wDir(idx(1)+1:end);
dirPath(strfind(dirPath, '\')) = '/';

IntegralS.wDirWsl = ['/mnt/d/',dirPath];
IntegralS.genParameterFile();
IntegralS.survivalExecutionCommand = ['wsl ',IntegralS.wDirWsl,'/', IntegralS.survivalParameterFileName, '.txt'];
IntegralS.execute();

%% ReadOut
IntegralS.readSingleIonLUT(1);
IntegralS.plotZDSingleIon(0);

%% Load LUT Matlab
fileNames = dir('D:\matRad_gitHubRemo_MKM_TOPAS\MKM\LUTs_rn35_RN_39_Beta005_500MeV_opt\*.mat');

for k=1:size(fileNames,1)
   LUT(k) = {load(strcat(fileNames(1).folder, filesep, fileNames(k).name))};
end

hold on;
plot(LUT{1}.E, LUT{1}.zD, '.-');
legend('Survival', 'Matlab');