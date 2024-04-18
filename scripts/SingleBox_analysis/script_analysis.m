%% load data
load('BOXPHANTOM_OAR_dis.mat');
matRad_cfg = MatRad_Config.instance();
saveDir = fullfile(matRad_cfg.matRadRoot, 'BOXPHANTOM_OAR_dist_analysis', '2mm');


load(fullfile(saveDir,'probQuantities.mat'), 'pln', 'cst', 'omega');
load(fullfile(saveDir, 'stf.mat'));

%% load weights


Case =7;
load(fullfile(saveDir, ['planResults_', num2str(Case), '.mat']), 'w', 'costFunctions');

%w = 12.1113*ones(size(w));
%% Analysis

RBE = 1.1;
scensToLoad = [1:10];%[1:newMultiScen.totNumScen];

[dij, newMultiScen, ~] = matRad_loadDijScenarios(ct,saveDir,'all', 'none');

physicalDose = cell(numel(scensToLoad),1);
physicalDose(:) = {zeros(dij.doseGrid.numOfVoxels,1)};

cstOnGrid = matRad_resizeCstToGrid(cst, dij.ctGrid.x, dij.ctGrid.y, dij.ctGrid.z, dij.doseGrid.x, dij.doseGrid.y, dij.doseGrid.z);

dvhDoseGrid = linspace(0,80,1000);
dvhs = [];
stdGrid = [];

stringLenght = 0;
scenCounter = 0;

for scenIdx=scensToLoad

    scenCounter = scenCounter +1;
    fprintf(repmat('\b',1,stringLenght));
    stringLenght = fprintf('\tComputing scenarios: %u/%u\n', scenCounter, numel(scensToLoad));
    
    [currDij, ~, ~] = matRad_loadDijScenarios(ct,saveDir,scenIdx,[],0);
    physicalDose{scenCounter} = currDij.physicalDose{1}*w*RBE*pln.numOfFractions;

end
    
% Get DVHs

dvhs = struct('name', [], 'volumePoints', []);
stringLenght = 0;
for scenCounter=scensToLoad
    fprintf(repmat('\b',1,stringLenght));
    stringLenght = fprintf('\t DVHs: %u/%u\n', scenCounter,numel(scensToLoad));

    tmpDVH = matRad_calcDVH(cstOnGrid,physicalDose{scenCounter},1,[],dvhDoseGrid);

    for k=1:numel(tmpDVH)
        dvhs(k).name = tmpDVH(k).name;
        dvhs(k).volumePoints = [dvhs(k).volumePoints; tmpDVH(k).volumePoints];

    end

end

% SDVHs

stringLenght = 0;

SDVH = struct('name', [], 'volumePoints', []);
if isempty(stdGrid)
    [tmpSDVH,tmpEx,tmpSD,stdGrid] = matRad_computeSTDVH(cstOnGrid,physicalDose,newMultiScen.scenWeight(scensToLoad)./sum(newMultiScen.scenWeight(scensToLoad)));
else
    [tmpSDVH,tmpEx,tmpSD,~] = matRad_computeSTDVH(cstOnGrid,physicalDose,newMultiScen.scenWeight(scensToLoad)./sum(newMultiScen.scenWeight(scensToLoad)),stdGrid);        
end


for k=1:numel(tmpSDVH)
    SDVH(k).name = tmpSDVH(k).name;
    SDVH(k).volumePoints = tmpSDVH(k).volumePoints;
end
SDVHdist = tmpSD;
SDVHdist_exp = tmpEx;

dvhExp = matRad_calcDVH(cstOnGrid,SDVHdist_exp,1,[],dvhDoseGrid);


for structIdx=1:size(cstOnGrid,1)
    allVoxels = arrayfun(@(scenStruct) scenStruct{1}, cstOnGrid{structIdx,4}, 'UniformOutput',false);
    structVoxels{structIdx} = unique([allVoxels{:}]);

    %this is equal to vTot/#voxels -> vTot/#Voxels = mean((SDVHdist(voxels)/pln.numOfFractions)^2)
    %meanSDVH{structIdx} = mean((SDVHdist(structVoxels{structIdx})/pln.numOfFractions).^2);
    meanSDVH{structIdx} = mean(SDVHdist(structVoxels{structIdx}));
end


%% Save


resultSaveDir = fullfile('BOXPHANTOM_OAR_dist_analysis', '2mm_results');
if ~exist(resultSaveDir, 'dir')
    mkdir(resultSaveDir);

end


fName = ['case_', num2str(Case), '_', num2str(numel(scensToLoad)), '.mat'];
save(fullfile(resultSaveDir, [fName]),'dvhExp', 'dvhs', 'SDVH', 'SDVHdist_exp', 'SDVHdist', 'dvhDoseGrid','stdGrid','meanSDVH', '-v7.3');
%% Vis


blue = [0.35 0.7 0.9];
orange = [0.9,0.6,0];
green = [0.4660 0.6740 0.1880];
red = [1 0 0];
colors = [blue; orange; green;red];

newCT = ct;
newCT.cubeDim = dij.doseGrid.dimensions;
newCT.cube{1} = interp3(dij.ctGrid.x, dij.ctGrid.y', dij.ctGrid.z,...
                        ct.cube{1}, ...
                        dij.doseGrid.x, dij.doseGrid.y', dij.doseGrid.z);

newCT.cubeHU{1} = interp3(dij.ctGrid.x, dij.ctGrid.y', dij.ctGrid.z,...
                        ct.cubeHU{1}, ...
                        dij.doseGrid.x, dij.doseGrid.y', dij.doseGrid.z);

newCT.x = dij.doseGrid.x;
newCT.y = dij.doseGrid.y;

newCT.z = dij.doseGrid.z;
figure;

movegui(gcf(), 'southwest');
leg = [];
structInclude = [2,3];
for structIdx =structInclude
    p(structIdx) = plot(gca(),dvhDoseGrid, dvhExp(structIdx).volumePoints, 'Color',colors(structIdx,:), 'Marker', 'none', 'LineStyle','-', 'LineWidth',1.5);
    hold on;

    for i=1:size(dvhs(structIdx).volumePoints,1)
         f = plot(gca(),dvhDoseGrid, dvhs(structIdx).volumePoints(i,:), 'Color',colors(structIdx,:), 'Marker', 'none', 'LineStyle','--', 'LineWidth',0.5);
         set(get(get(f,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end

    leg = [leg, {dvhs(structIdx).name}];
end
grid on;

legend(leg);

%% Dose dist
doseCube = reshape(physicalDose{3}, dij.doseGrid.dimensions);

%doseCube = reshape(SDVHdist_exp, dij.doseGrid.dimensions);

slice = 121;
profile = 121;
figure;
movegui(gcf(),'southwest');
tiledlayout(1,2);
nexttile;
matRad_plotSliceWrapper(gca(),newCT,cstOnGrid,1,doseCube,3,slice);
xline(profile);

nexttile;

for scenIdx=1:numel(scensToLoad)
    doseCube = reshape(physicalDose{scenIdx}, dij.doseGrid.dimensions);
    hold on;
    plot(dij.doseGrid.y, squeeze(doseCube(:,profile,slice)), '.-');
end
%% Cost functions

iters = 1:numel(costFunctions(1).values);
figure;
movegui(gcf(), 'southwest');
leg = [];
for costIdx=1:numel(costFunctions)
    semilogy(iters,costFunctions(costIdx).values, '.-');
    hold on;
    leg = [leg, {costFunctions(costIdx).name}];
end

grid on;
legend(leg);

%% STD distribution
stdCube = reshape(SDVHdist, dij.doseGrid.dimensions);

figure;
movegui(gcf(),'southwest');
tiledlayout(1,2);
nexttile;
matRad_plotSliceWrapper(gca(),newCT,cstOnGrid,1,stdCube,3,slice,[],[],[],[],[0,30]);
xline(profile);

nexttile;

for scenIdx=1:1
    %doseCube = reshape(physicalDose{scenIdx}, dij.doseGrid.dimensions);
    hold on;
    plot(dij.doseGrid.y, squeeze(stdCube(:,profile,slice)), '.-');
end


%% Plot SDVHs
figure;
movegui(gcf(), 'southwest');
leg = [];
structInclude = [2,3];
for structIdx =structInclude
    p(structIdx) = plot(gca(),stdGrid, SDVH(structIdx).volumePoints, 'Color',colors(structIdx,:), 'Marker', 'none', 'LineStyle','-', 'LineWidth',1.5);
    hold on;
    leg = [leg, {SDVH(structIdx).name}];

    allVoxels = arrayfun(@(scenStruct) scenStruct{1}, cstOnGrid{structIdx,4}, 'UniformOutput',false);
    structVoxels{structIdx} = unique([allVoxels{:}]);

end
grid on;

legend(leg);

% 
% varStruct = 2;
% dOmega = w'*omega{varStruct}*RBE*RBE;
% vTot = dOmega*w;
% 
% 
% meanVTot = (vTot/numel(structVoxels{varStruct}));
% varParameter = 0.02;
% 
% figure;
% xline(meanVTot);
% xline(varParameter, 'Color','r');
% 
% xlim([0,max([varParameter, meanVTot])+1]);
% 
% %figure;
% meanSDVH = mean((SDVHdist(structVoxels{2})/pln.numOfFractions).^2);
% 
% xline(meanSDVH);





%% Compare results
casesToInclude = [1,5];
variables = {'dvhExp','dvhs', 'SDVH', 'SDVHdist_exp', 'SDVHdist','dvhDoseGrid', 'stdGrid'};
clear case_1 case_2 case_3 case_4;
clear(variables{:});
for caseIdx=casesToInclude

    load(fullfile(resultSaveDir, ['case_', num2str(caseIdx), '_',num2str(numel(scensToLoad)), '.mat']));

    for varIdx=1:numel(variables)
        currCase.(variables{varIdx}) = eval(variables{varIdx});
       
    end
    eval(['case_', num2str(caseIdx),' = currCase;']);

    clear currCase;
end

%% compare DVHs
fontSize = 17;

figure;
movegui(gcf(), 'southwest');
currPosition = get(gcf(), 'Position');
figureWidth = 700;
aspectRatio = 3/4;
figureHeight = figureWidth/aspectRatio;

set(gcf(), 'Position', [currPosition(1), currPosition(2), figureWidth, figureHeight]);


%tiledlayout(2,1);
structCounter=0;
structInclude = [2,3];
for structIdx = structInclude
    %nexttile;
    structCounter = structCounter+1;
    subplot(2,1,structCounter);
    ax = gca();
    currAxPos = get(ax, 'Position');
    %Position is normalized
    axPosition(1) = 0.09;  %H position
    axPosition(2) = 0.055 + (2-structCounter)*0.5; %V position
    axWidth = 0.85; % H dimension
    axHeight = 0.40; % V dimension
    set(ax, 'Position', [axPosition, axWidth, axHeight]);

    leg = [];
    caseCounter = 0;
    for caseIdx=casesToInclude
        caseCounter = caseCounter +1;
        currCase = eval(['case_', num2str(caseIdx)]);
        p(structIdx) = plot(gca(),currCase.dvhDoseGrid, currCase.dvhExp(structIdx).volumePoints, 'Color',colors(caseCounter,:), 'Marker', 'none', 'LineStyle','-', 'LineWidth',1.5);
        hold on;
    
        for i=1:size(currCase.dvhs(structIdx).volumePoints,1)
            f = plot(gca(),currCase.dvhDoseGrid, currCase.dvhs(structIdx).volumePoints(i,:), 'Color',colors(caseCounter,:), 'Marker', 'none', 'LineStyle','--', 'LineWidth',0.5);
             set(get(get(f,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
    
        leg = [leg, {['case_', num2str(caseIdx)]}];
        
    end
    grid on;
    legend(leg);
end


ax = get(gcf(), 'Children');

% Child #4 is the top subplot
cNumb = 4;
ax(cNumb).Title.String = 'Target DVH';
ax(cNumb).Title.FontSize = fontSize;
ax(cNumb).YLabel.String = 'Volume [%]';
ax(cNumb).YLabel.FontSize = fontSize;
ax(cNumb).XLabel.String = 'Dose [Gy]';
ax(cNumb).XLabel.FontSize = fontSize;
ax(cNumb).XLim = [50, 70];
% Child #2 is the bottom subplot
cNumb = 2;
ax(cNumb).Title.String = 'OAR DVH';
ax(cNumb).Title.FontSize = fontSize;
ax(cNumb).YLabel.String = 'Volume [%]';
ax(cNumb).YLabel.FontSize = fontSize;
ax(cNumb).XLabel.String = 'Dose [Gy]';
ax(cNumb).XLabel.FontSize = fontSize;


% Legend plot top
cNumb = 3;
ax(cNumb).String = {'scenario-free', 'probabilistic'};
ax(cNumb).FontSize = fontSize;

% Legend bottom
cNumb = 1;
ax(cNumb).String = {'scenario-free', 'probabilistic'};
ax(cNumb).FontSize = fontSize;



%% SDVH comparison
% compare DVHs
figure;
movegui(gcf(), 'southwest');
currPosition = get(gcf(), 'Position');
figureWidth = 700;
aspectRatio = 3/4;
figureHeight = figureWidth/aspectRatio;

set(gcf(), 'Position', [currPosition(1), currPosition(2), figureWidth, figureHeight]);

%tiledlayout(2,1);

structInclude = [2,3];
structCounter = 0;
for structIdx = structInclude
    structCounter = structCounter+1;
    %nexttile;
    subplot(2,1,structCounter);
    ax = gca();
    currAxPos = get(ax, 'Position');
    %Position is normalized
    axPosition(1) = 0.09;  %H position
    axPosition(2) = 0.05 + (2-structCounter)*0.5; %V position
    axWidth = 0.85; % H dimension
    axHeight = 0.40; % V dimension
    set(ax, 'Position', [axPosition, axWidth, axHeight]);

    leg = [];
    caseCounter = 0;
    for caseIdx=casesToInclude
        caseCounter = caseCounter +1;
        currCase = eval(['case_', num2str(caseIdx)]);
        p(structIdx) = plot(gca(),currCase.stdGrid, currCase.SDVH(structIdx).volumePoints, 'Color',colors(caseCounter,:), 'Marker', 'none', 'LineStyle','-', 'LineWidth',1.5);
        hold on;
    
        leg = [leg, {['case_', num2str(caseIdx)]}];
        
    end
    grid on;
    legend(leg);
end

%nexttile(1);
% Get axes
ax = get(gcf(), 'Children');

% Child #4 is the top subplot
cNumb = 4;
ax(cNumb).Title.String = 'Target SDVH';
ax(cNumb).Title.FontSize = fontSize;
ax(cNumb).YLabel.String = 'Volume [%]';
ax(cNumb).YLabel.FontSize = fontSize;
ax(cNumb).XLabel.String = 'Std [Gy]';
ax(cNumb).XLabel.FontSize = fontSize;
ax(cNumb).XLim = [0, 17];
% Child #2 is the bottom subplot
cNumb = 2;
ax(cNumb).Title.String = 'OAR SDVH';
ax(cNumb).Title.FontSize = fontSize;
ax(cNumb).YLabel.String = 'Volume [%]';
ax(cNumb).YLabel.FontSize = fontSize;
ax(cNumb).XLabel.String = 'Std [Gy]';
ax(cNumb).XLabel.FontSize = fontSize;


% Legend plot top
cNumb = 3;
ax(cNumb).String = {'scenario-free', 'probabilistic'};
ax(cNumb).FontSize = fontSize;

% Legend bottom
cNumb = 1;
ax(cNumb).String = {'scenario-free', 'probabilistic'};
ax(cNumb).FontSize = fontSize;

