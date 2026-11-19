matRad_cfg = MatRad_Config.instance(); 
outDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','doseCalculation','analyticalEngines');
outFile = fullfile(outDir,'protons_insert_analytical.mat');

matRad_cfg = MatRad_Config.instance();
load(fullfile(outFile));

%% PLOT TO TICKZ
% saveDir = fullfile(matRad_cfg.primaryUserFolder, 'dataPaper', 'Figures', 'doseCalculation', 'analyticalEngines', 'protons_insert', 'Hong');
% saveDir = fullfile(matRad_cfg.primaryUserFolder, 'dataPaper', 'Figures', 'doseCalculation', 'analyticalEngines', 'protons_insert', 'fineSampling');
saveDir = fullfile(matRad_cfg.primaryUserFolder, 'dataPaper', 'Figures', 'doseCalculation', 'analyticalEngines', 'protons_insert', 'Bortfeld');

cstPlot = cst;
% doseCube = doseHong;
% doseCube = doseFineSampling;
doseCube = doseBortfeld;

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
                  'dose', double(doseCube), ...
                  'plane', 3, ...
                  'slice',27,...
                  'doseIsoLevels', 1.8e-3*[50,70,80,95,105]./100,...
                  'ctWindow', [lowCTLim,higCtLim]+level, ...
                  'alpha', 0.55,...
                  'doseWindow', [0,1.8e-3],...
                  'voiSelection', [1 1 0]);...
                  % 'doseColorMap', diffMap);

% xlim([75,155]-5);
% ylim([75,155]-5);

title('');

if savePlots
    fileName = sprintf('DoseDist.tex');
    
    figSaveDir = fullfile(saveDir);
    
    if ~exist(figSaveDir, 'dir')
        mkdir(figSaveDir);
    end
    
    matRad_saveFigToTikz(f, figSaveDir, fileName, hContour, hIsoDose);
end