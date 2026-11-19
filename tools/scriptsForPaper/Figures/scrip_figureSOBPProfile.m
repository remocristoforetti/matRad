clear all;
matRad_cfg = MatRad_Config.instance();

% outDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','doseCalculation','analyticalEngines');
% outFile = fullfile(outDir,'protons_SOBP_analytical_2.mat');

% outDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','doseCalculation','analyticalEngines');
% outFile = fullfile(outDir,'protons_SOBPslabs_analytical_2.mat');
% outFile = fullfile(outDir,'carbon_SOBPslabs_analytical_2.mat');
% outFile = fullfile(outDir,'helium_SOBPslabs_analytical_2.mat');
% outFile = fullfile(outDir,'oxygen_SOBPslabs_analytical_2.mat');

outDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','doseCalculation','MCengines');
outFile = fullfile(outDir,'protons_SOBPslabs_MC_2.mat');
load(fullfile(outFile));

%% PLOT TO TICKZ
% saveDir = fullfile(matRad_cfg.primaryUserFolder, 'dataPaper', 'Figures', 'doseCalculation', 'analyticalEngines', 'slab_SOBP_2');
saveDir = fullfile(matRad_cfg.primaryUserFolder, 'dataPaper', 'Figures', 'doseCalculation', 'MCengines', 'slab_SOBP_2');

particle = 'protons';
% particle = 'carbon';
% particle = 'helium';
% particle = 'oxygen';

% fileName = sprintf('%s_HONG.txt', particle);
fileName = sprintf('%s_slabMC_TOPAS.txt', particle);
% fileName = sprintf('%s_slabMC_MCsquare.txt', particle);
% fileName = sprintf('%s_slabMC_FRED.txt', particle);

% doseCube = doseSOBP;
doseCube = doseTopasIon;
% doseCube = doseMCsquare;
% doseCube = doseFred;

savePlots = true;

% average over a few voxels around the central beam axis

isoIdx = [25,1,25];%matRad_world2cubeIndex(pln.propStf.isoCenter,ct); % [yIdx xIdx zIdx]

yIdx = isoIdx(1) + [-3:3];
zIdx = isoIdx(3) + [-3:3];

profile = squeeze(mean(doseCube(yIdx,:,zIdx), [1,3]));

f = figure;
hProfile = plot(ct.x- ct.x(1), profile, 'LineWidth', 2);
grid on;
xlabel('x [mm]');
ylabel('physical dose [Gy]');
title('');

if savePlots

    figSaveDir = fullfile(saveDir);

    if ~exist(figSaveDir, 'dir')
        mkdir(figSaveDir);
    end

    opts.append     = false;
    opts.forgetPlot = false;
    matRad_linesToTikzFile(hProfile, fullfile(figSaveDir, fileName), opts, 'profileLineWidth');
end