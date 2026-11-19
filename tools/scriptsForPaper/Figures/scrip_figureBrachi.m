clear all;
matRad_cfg = MatRad_Config.instance();
load('PROSTATE.mat');

outDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','brachytherapy');

outFile = fullfile(outDir,'brachy_prostate.mat');
load(fullfile(outFile));


%% PLOT TO TICKZ
saveDir = fullfile(matRad_cfg.primaryUserFolder, 'dataPaper', 'Figures', 'brachitherapy', 'DoseDist');
cstPlot = cst;

doseCube = doseBrachy;

savePlots = true;

lowCTLim =-400;
higCtLim =1140;

level = -580;

f = figure;
movegui(gcf(), 'southwest');
f.Position(1:2) = f.Position(1:2);%-800;
f.Position(3:4) = [900,900];
[hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSlice(ct, 'cst', cstPlot, ...
                 'axesHandle', gca(),...
                  'dose', double(doseCube), ...
                  'plane', 3, ...
                  'slice',33,...
                  'doseIsoLevels', 1.8e-3*[50,70,80,95,105]./100,...
                  'ctWindow', [lowCTLim,higCtLim]+level, ...
                  'alpha', 0.55,...
                  'doseWindow', [0,15],...
                  'voiSelection', [1 0 0 0 0 1 0 1 1 0]);...
                  % 'doseColorMap', diffMap);

xlim([55,125]-3);
ylim([55,125]);

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
saveDir = fullfile(matRad_cfg.primaryUserFolder, 'dataPaper', 'Figures', 'brachitherapy');

dvh = matRad_calcDVH(cst, doseBrachy,[],linspace(0,20,100));

structIdx = [1,4,6];
structName = {'Rectum','Prostate', 'PTV'};

for i=1:numel(structIdx)
    
    cdvh = dvh(structIdx(i));

    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end

    fileName = fullfile(saveDir, sprintf('dvh_%s.txt', structName{i}));

    matRad_writeTickzDVH(cdvh, fileName);
end