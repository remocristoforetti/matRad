clear all;
matRad_cfg = MatRad_Config.instance();
load(fullfile(matRad_cfg.primaryUserFolder,'patients','TG166_brain.mat'),'ct','cst');
outDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','bioModels');


% outFile = fullfile(outDir,'protons_TG166brain_TOPAS_MCN.mat');
% outFile = fullfile(outDir,'protons_TG166brain_MCN.mat');
% load(fullfile(outFile));

outFileHong = fullfile(outDir,'protons_TG166brain_MCN.mat');
outFileTopas = fullfile(outDir,'protons_TG166brain_TOPAS_MCN.mat');

load(fullfile(outFileHong));
load(fullfile(outFileTopas));

%% PLOT TO TICKZ
% saveDir = fullfile(matRad_cfg.primaryUserFolder, 'dataPaper', 'Figures', 'bioModels', 'TG166Brain','TOPAS', 'DoseDist');
saveDir = fullfile(matRad_cfg.primaryUserFolder, 'dataPaper', 'Figures', 'bioModels', 'TG166Brain','Gamma', 'DoseDist');
cstPlot = cst;

doseCube =matRad_gammaIndex(RBExDose, RBExDoseTopas, [ct.resolution.x, ct.resolution.y, ct.resolution.z], [3,3]);    ;

% doseCube = interp3(doseGrid.x, doseGrid.y', doseGrid.z, doseCube, ct.x, ct.y',ct.z);

savePlots = false;

lowCTLim =-400;
higCtLim =1140;

level = -580;

f = figure;
movegui(gcf(), 'southwest');
f.Position(1:2) = f.Position(1:2);%-800;
f.Position(3:4) = [900,900];
[hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSlice(ct, 'cst', cstPlot, ...
                 'axesHandle', gca(),...
                  'dose', doseCube, ...
                  'plane', 3, ...
                  'slice',37,...
                  'doseIsoLevels', 60*[50,70,80,95,105]./100,...
                  'ctWindow', [lowCTLim,higCtLim]+level, ...
                  'alpha', 0.55,...
                  ...% 'doseWindow', [0,70],...
                  'doseWindow', [0,2],...
                  'voiSelection', [0 0 1 1 1 0 0 0 0 0 0 0 0 1 0 0],...
                  'doseColorMap', gammaIndex);

xlim([150,350]+3);
ylim([150,350]+10);

title('');

if savePlots
    fileName = sprintf('DoseDist.tex');
    
    figSaveDir = fullfile(saveDir);
    
    if ~exist(figSaveDir, 'dir')
        mkdir(figSaveDir);
    end
    
    hContour(cellfun(@isempty,hContour)) = [];
    matRad_saveFigToTikz(f, figSaveDir, fileName, hContour, hIsoDose);
end

%% DVH
downSampling = 3;

saveDir = fullfile(matRad_cfg.primaryUserFolder, 'dataPaper', 'Figures','robustness', 'P34');
cData = analysis.(plan);

dvhIdxs = [1:downSampling:1000];

structIdx = [1,2,3];
structName = {'GTV','Lung', 'PTV'};
for i=1:numel(structIdx)
    
    cdvh = cData.dvhNominal(structIdx(i));

    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end

    fileName = fullfile(saveDir, sprintf('%s_dvh_%s_nominal.txt',plan, structName{i}));

    matRad_writeTickzDVH(cdvh, fileName);

    dvh_l = cData.dvhPercentile_5(structIdx(i));
    dvh_l.doseGrid     = dvh_l.doseGrid(dvhIdxs);
    dvh_l.volumePoints = dvh_l.volumePoints(dvhIdxs);
    
    dvh_h = cData.dvhPercentile_95(structIdx(i));
    
    dvh_h.doseGrid     = dvh_h.doseGrid(dvhIdxs);
    dvh_h.volumePoints = dvh_h.volumePoints(dvhIdxs);
    
    fileName = fullfile(saveDir, sprintf('%s_dvh_%s_band.txt',plan, structName{i}));
    matRad_writeTickzDVH_Fill(dvh_h, dvh_l, fileName);
end