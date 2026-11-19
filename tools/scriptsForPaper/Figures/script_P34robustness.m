clear all;
matRad_cfg = MatRad_Config.instance();
load('Lung_Patient_34.mat','ct','cst');
outDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','robustness', 'P34');


outFile = fullfile(outDir,'robustnessAnalysis.mat');
load(fullfile(outFile));

plan = 'nominal';
%% PLOT TO TICKZ
saveDir = fullfile(matRad_cfg.primaryUserFolder, 'dataPaper', 'Figures', 'robustness', 'P34', 'DoseDist');
cstPlot = cst;

doseCube = analysis.(plan).doseNominal;

doseCube = interp3(doseGrid.x, doseGrid.y', doseGrid.z, doseCube, ct.x, ct.y',ct.z);

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
                  'slice',60,...
                  'doseIsoLevels', 60*[50,70,80,95,105]./100,...
                  'ctWindow', [lowCTLim,higCtLim]+level, ...
                  'alpha', 0.55,...
                  'doseWindow', [0,70],...
                  'voiSelection', [1 1 1 0 0 0 0 0 0 0]);...
                  % 'doseColorMap', diffMap);

xlim([300,500]-50);
ylim([300,500]-200);

title('');

if savePlots
    fileName = sprintf('DoseDist.tex');
    
    figSaveDir = fullfile(saveDir);
    
    if ~exist(figSaveDir, 'dir')
        mkdir(figSaveDir);
    end
    
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