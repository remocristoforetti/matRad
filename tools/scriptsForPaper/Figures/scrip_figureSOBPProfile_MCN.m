clear all;
matRad_cfg = MatRad_Config.instance();

% outDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','bioModels');
% outFile = fullfile(outDir,'protons_SOBP_bioModels_analytical_2.mat');

outDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','bioModels');
outFile = fullfile(outDir,'protons_SOBP_bioModels_TOPAS_MCN.mat');
% outFile = fullfile(outDir,'carbon_SOBPslabs_analytical_2.mat');
% outFile = fullfile(outDir,'helium_SOBPslabs_analytical_2.mat');
% outFile = fullfile(outDir,'oxygen_SOBPslabs_analytical_2.mat');

load(fullfile(outFile));

%% PLOT TO TICKZ
saveDir = fullfile(matRad_cfg.primaryUserFolder, 'dataPaper', 'Figures', 'doseCalculation', 'bioModels');
fileName = sprintf('LET_TOPAS.txt');

% doseCube = RBExDose.LSM;
% doseCube = physicalDoseTopas;
doseCube = LETTopas;

savePlots = false;

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