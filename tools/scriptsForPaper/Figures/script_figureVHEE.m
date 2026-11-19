clear all;
matRad_cfg = MatRad_Config.instance();
load('LIVER.mat');

outDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','VHEE');

outFile = fullfile(outDir,'VHEE_liver.mat');
load(fullfile(outFile));


%% PLOT TO TICKZ
saveDir = fullfile(matRad_cfg.primaryUserFolder, 'dataPaper', 'Figures', 'VHEE', 'DoseDist');
cstPlot = cst;

doseCube = doseVHEE*30;

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
                  'slice',130,...
                  'doseIsoLevels', 50*[50,70,80,95,105]./100,...
                  'ctWindow', [lowCTLim,higCtLim]+level, ...
                  'alpha', 0.55,...
                  'doseWindow', [0,50],...
                  'voiSelection', [0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0]);...
                  % 'doseColorMap', diffMap);

xlim([50,130]);
ylim([50,130]);

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
saveDir = fullfile(matRad_cfg.primaryUserFolder, 'dataPaper', 'Figures', 'VHEE');

dvh = matRad_calcDVH(cst, doseVHEE*30,[],linspace(0,50,300));

structIdx = [10,11,16];
structName = {'Heart','SpinalChord', 'PTV'};

for i=1:numel(structIdx)
    
    cdvh = dvh(structIdx(i));

    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end

    fileName = fullfile(saveDir, sprintf('dvh_%s.txt', structName{i}));

    matRad_writeTickzDVH(cdvh, fileName);
end