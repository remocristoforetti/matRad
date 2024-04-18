%% load data
load('4D_BOXPHANTOM_OAR_dist_heavy.mat');

% cubeIdx = [1,3,5];
% cubes = [];
% cubesHU = [];
% 
% for i=1:3
% 
%     cubes{i}   = ct.cube{cubeIdx(i)};
%     cubesHU{i} = ct.cubeHU{cubeIdx(i)};
% end
% 
% ct.cube = cubes;
% ct.cubeHU = cubesHU;
% 
% ct.numOfCtScen = 3;


matRad_cfg = MatRad_Config.instance();
saveDir = fullfile(matRad_cfg.matRadRoot, '4D_BOXPHANTOM_OAR_dist_heavy_analysis', '2mm');


load(fullfile(saveDir,'probQuantities.mat'), 'pln', 'cst');
load(fullfile(saveDir, 'stf.mat'));

%% load weights


Case =4;
load(fullfile(saveDir, ['planResults_', num2str(Case), '.mat']),'prob_info', 'w', 'costFunctions');
%load(fullfile(saveDir, ['planResults_phase_', num2str(Case), '.mat']), 'w', 'costFunctions');
%load(fullfile(matRad_cfg.matRadRoot, '4D_BOXPHANTOM_OAR_dist_heavy_analysis', 'ITV', ['plan_results_', num2str(Case), '.mat']), 'prob_info', 'w', 'costFunctions');

%% Analysis

RBE = 1.1;
scensToLoad = [1:30];%[1:newMultiScen.totNumScen];

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
    
    [currCtIdx, ~,~] = ind2sub(size(dij.physicalDose),scenIdx);
    [currDij, ~, ~] = matRad_loadDijScenarios(ct,saveDir,scenIdx,[],0);
    physicalDose{scenCounter} = currDij.physicalDose{currCtIdx}*w*RBE*pln.numOfFractions;

end
    
% Get DVHs

dvhs = struct('name', [], 'volumePoints', []);
stringLenght = 0;
for scenCounter=scensToLoad
    fprintf(repmat('\b',1,stringLenght));
    stringLenght = fprintf('\t DVHs: %u/%u\n', scenCounter,numel(scensToLoad));


    [currCtIdx, ~,~] = ind2sub(size(dij.physicalDose),scenCounter);
    tmpDVH = matRad_calcDVH(cstOnGrid,physicalDose{scenCounter},currCtIdx,[],dvhDoseGrid);

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
dvhExp = matRad_calcDVH(cstOnGrid,SDVHdist_exp,[1:ct.numOfCtScen],[],dvhDoseGrid);

structVoxels = {};
for structIdx=1:size(cstOnGrid,1)
    allVoxels = arrayfun(@(scenStruct) scenStruct{1}', cstOnGrid{structIdx,4}, 'UniformOutput',false);
    structVoxels{structIdx} = unique([allVoxels{:}])';

    %this is equal to vTot/#voxels -> vTot/#Voxels = mean((SDVHdist(voxels)/pln.numOfFractions)^2)
    %meanSDVH{structIdx} = mean((SDVHdist(structVoxels{structIdx})/pln.numOfFractions).^2);

    meanSDVH{structIdx} = mean(SDVHdist(structVoxels{structIdx}));
end

%% Save

resultSaveDir = fullfile('4D_BOXPHANTOM_OAR_dist_heavy_analysis', '2mm_results');
if ~exist(resultSaveDir, 'dir')
    mkdir(resultSaveDir);


end

fName = ['case_', num2str(Case), '_', num2str(numel(scensToLoad)), '.mat'];
%fName = ['phase_case_', num2str(Case), '_', num2str(numel(scensToLoad)), '.mat'];

save(fullfile(resultSaveDir, [fName]),'dvhExp', 'dvhs', 'SDVH', 'SDVHdist_exp', 'SDVHdist', 'dvhDoseGrid','stdGrid','meanSDVH', '-v7.3');
%% Vis



blue = [0.35 0.7 0.9];


orange = [0.9,0.6,0];
green = [0.4660 0.6740 0.1880];
red = [1 0 0];
colors = [blue; orange; green;red];

newCT = ct;
newCT.cubeDim = dij.doseGrid.dimensions;
newCT.cube = cellfun(@(cube) interp3(dij.ctGrid.x, dij.ctGrid.y', dij.ctGrid.z,...
                        cube, ...
                        dij.doseGrid.x, dij.doseGrid.y', dij.doseGrid.z), ct.cube, 'UniformOutput',false);

newCT.cubeHU = cellfun(@(cube) interp3(dij.ctGrid.x, dij.ctGrid.y', dij.ctGrid.z,...
                        cube, ...
                        dij.doseGrid.x, dij.doseGrid.y', dij.doseGrid.z), ct.cubeHU, 'UniformOutput',false);

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
doseCube = reshape(SDVHdist_exp, dij.doseGrid.dimensions);
%doseCube = reshape(physicalDose{3}, dij.doseGrid.dimensions);

%doseCube = reshape(SDVHdist_exp, dij.doseGrid.dimensions);
cubeIdx = 3;
slice = 121;
profile = 121;
figure;
movegui(gcf(),'southwest');
tiledlayout(1,3);
nexttile;
matRad_plotSliceWrapper(gca(),newCT,cstOnGrid,cubeIdx,doseCube,3,slice);
xline(profile);

nexttile;

for scenIdx=1:numel(scensToLoad)
    doseCube = reshape(physicalDose{scenIdx}, dij.doseGrid.dimensions);

    hold on;
    plot(dij.doseGrid.y, squeeze(doseCube(:,profile,slice)), '.-');
end

nexttile;
for scenIdx=1:numel(scensToLoad)

    doseCube = reshape(physicalDose{scenIdx}, dij.doseGrid.dimensions);
    hold on;
    plot(dij.doseGrid.y, squeeze(doseCube(profile,:,slice)), '.-');
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
end

grid on;

legend(leg);



%% Compare results
casesToInclude = [1];
variables = {'dvhExp','dvhs', 'SDVH', 'SDVHdist_exp', 'SDVHdist', 'dvhDoseGrid', 'stdGrid'};
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

for caseIdx=casesToInclude

    load(fullfile(resultSaveDir, ['phase_case_', num2str(caseIdx), '_',num2str(numel(scensToLoad)), '.mat']));

    for varIdx=1:numel(variables)
        currCase.(variables{varIdx}) = eval(variables{varIdx});
       
    end
    eval(['phase_case_', num2str(caseIdx),' = currCase;']);

    clear currCase;
end



%% compare DVHs phase/all
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

    for caseIdx=casesToInclude
        caseCounter = caseCounter +1;
        currCase = eval(['phase_case_', num2str(caseIdx)]);
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
ax(cNumb).Title.String = ['Target DVH Case ', num2str(casesToInclude)];
ax(cNumb).Title.FontSize = fontSize;
ax(cNumb).YLabel.String = 'Volume [%]';
ax(cNumb).YLabel.FontSize = fontSize;
ax(cNumb).XLabel.String = 'Dose [Gy]';
ax(cNumb).XLabel.FontSize = fontSize;
ax(cNumb).XLim = [40, 70];
% Child #2 is the bottom subplot
cNumb = 2;
ax(cNumb).Title.String = ['OAR DVH Case ', num2str(casesToInclude)];
ax(cNumb).Title.FontSize = fontSize;
ax(cNumb).YLabel.String = 'Volume [%]';
ax(cNumb).YLabel.FontSize = fontSize;
ax(cNumb).XLabel.String = 'Dose [Gy]';
ax(cNumb).XLabel.FontSize = fontSize;


% Legend plot top
cNumb = 3;
ax(cNumb).String = {'all', 'phase'};
ax(cNumb).FontSize = fontSize;

% Legend bottom
cNumb = 1;
ax(cNumb).String = {'all', 'phase'};
ax(cNumb).FontSize = fontSize;


%% SDVH comparison
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


    for caseIdx=casesToInclude
        caseCounter = caseCounter +1;
        currCase = eval(['phase_case_', num2str(caseIdx)]);
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
ax(cNumb).XLim = [0, 2];
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
ax(cNumb).FontSize = fontSize;%% Vis




%% gif
figure;
%movegui(gcf(), 'southwest');
for scenIdx=scensToLoad
    clf;
    %imagesc(ct.cubeHU{i}(:,:,80));
    %matRad_plotVoiContourSlice(gca(),cst,ct,i,[],3,80);
    [currCtIdx, ~,~] = ind2sub(size(dij.physicalDose),scenIdx);

    doseCube = reshape(physicalDose{scenIdx}, dij.doseGrid.dimensions);
    matRad_plotSliceWrapper(gca(),newCT,cstOnGrid,currCtIdx,doseCube,3,121);
    pause(0.5);
    fram(scenIdx) = getframe(gcf);
    im{scenIdx} = frame2im(fram(scenIdx));
end

dirName = fullfile(saveDir, 'gifs');


if ~exist(dirName, 'dir')
    mkdir(dirName);
end

fName = fullfile(dirName, ['testGIF', '.gif']);

for scenIdx=scensToLoad
    [img,map] = rgb2ind(im{scenIdx},256);

    if scenIdx==1
        imwrite(img, map, fName, 'gif',"LoopCount",Inf,"DelayTime",0.5);
    else
        imwrite(img,map,fName,"gif","WriteMode","append","DelayTime",0.5);
    end

end

% %% %% Box dose dist
% load('4D_BOXPHANTOM_OAR_dist_heavy.mat');
% 
% newCT = ct;
% useDoseGrid = 1;
% if useDoseGrid
%     newCT.cubeDim = dij.doseGrid.dimensions;
%     newCT.cubeHU = cellfun(@(cube) interp3(dij.ctGrid.x, dij.ctGrid.y', dij.ctGrid.z,...
%                             cube, ...
%                             dij.doseGrid.x, dij.doseGrid.y', dij.doseGrid.z), ct.cubeHU, 'UniformOutput',false);
% 
%     newCT.x = dij.doseGrid.x;
%     newCT.y = dij.doseGrid.y;
%     newCT.z = dij.doseGrid.z;
% 
%     cstOnGrid = matRad_resizeCstToGrid(ITVcst,dij.ctGrid.x, dij.ctGrid.y, dij.ctGrid.z,newCT.x, newCT.y, newCT.z);
% else
%     cstOnGrid = cst;
% end


%% Plot
figureWidth = 900;
aspectRatio = 4/3;
figureHeight = figureWidth/aspectRatio;

fontSize = 17;
slice = 121;

figure;
movegui(gcf(),'southwest');
currPosition = get(gcf(), 'Position');
set(gcf(), 'Position', [currPosition(1), currPosition(2), figureWidth, figureHeight]);

%tiledlayout(numel(cases), 2);

ctIdx=1;
    subplot(2,2,1);
    cube = reshape(physicalDose{ctIdx}, newCT.cubeDim);

    matRad_plotSliceWrapper(gca(),newCT, cstOnGrid, ctIdx, cube,3, slice, [], [],[],[],[0,70], [], [0 1 1 1 0 0 0 1 0 0 0 0 1 0 0]);
    title('Dose CT scenario 1', 'FontSize', fontSize+3, 'FontWeight','bold');
    ylabel('');
 
ctIdx = 6;
    subplot(2,2,2);
    cube = reshape(physicalDose{ctIdx}, newCT.cubeDim);
    matRad_plotSliceWrapper(gca(),newCT, cstOnGrid, ctIdx, cube,3, slice, [], [],[],[],[0,70], [], [0 1 1 1 0 0 0 0 0 0 0 0 1 0 1]);
    title('Dose CT scenario 5', 'FontSize',fontSize+3);

ctIdx = 10;
    currCase = eval([casePreName{caseIdx}, num2str(cases(caseIdx))]);
    %nexttile;
    subplot(2,2,3);
    cube = reshape(physicalDose{ctIdx}, newCT.cubeDim);
    matRad_plotSliceWrapper(gca(),newCT, cstOnGrid, ctIdx, cube,3, slice, [], [],[],[],[0,70], [], [0 1 1 1 0 0 0 1 0 0 0 0 1 0 0]);
    title('Dose CT scenario 10', 'FontSize', fontSize+3, 'FontWeight','bold');
    ylabel('');

    subplot(2,2,4);
    cube = reshape(SDVHdist, newCT.cubeDim);
    matRad_plotSliceWrapper(gca(),newCT, cstOnGrid, 6, cube,3, slice, [], [],[],[],[0,30], [], [0 1 1 1 0 0 0 0 0 0 0 0 1 0 1]);
    title('SD Distirbution', 'FontSize', fontSize+3, 'FontWeight','bold');


    ax = get(gcf(), 'Children');
% Child #4 is the top subplot
cNumb = [2,4,6,8];
for cNumbIdx=cNumb
    ax(cNumbIdx).XTickLabel = '';
    ax(cNumbIdx).XLabel.String = '';

    ax(cNumbIdx).YTickLabel = '';
end

ax(2).YLabel.String = '';
ax(6).YLabel.String = '';


%Lower right
ax(2).Position = [0.45, 0, 0.5, 0.5];
ax(2).OuterPosition = ax(2).Position;

% Lower left
ax(4).Position = [0, 0, 0.5, 0.5];
ax(4).OuterPosition = ax(4).Position;

% Upper right
ax(6).Position = [0.45, 0.47, 0.5, 0.53];
ax(6).OuterPosition = ax(6).Position;

%Upper left
ax(8).Position = [0, 0.47, 0.5, 0.53];
ax(8).OuterPosition = ax(8).Position;


% Colorbars
%Lower right
colorbarsIndexes = [1,3,5,7];
for cIdx=colorbarsIndexes
    ax(cIdx).Label.String = '[Gy]';
    ax(cIdx).Label.FontSize = 17;
    ax(cIdx).FontSize = 12;
end