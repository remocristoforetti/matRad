%% load data
load('C:\r408i_data\r408i_data\CTDatasetMotion\102_HM10395_333_newCst.mat');
matRad_cfg = MatRad_Config.instance();
saveDir = fullfile(matRad_cfg.matRadRoot, 'Patient_1_analysis', '2mm_30_patient');


load(fullfile(saveDir, 'stf.mat'));

%% load weights


Case =101;
load(fullfile(saveDir, ['planResults_', num2str(Case), '.mat']), 'w', 'costFunctions', 'cst', 'pln');


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


resultSaveDir = fullfile(matRad_cfg.matRadRoot, 'Patient_1_analysis', '2mm_30_patient_results');
if ~exist(resultSaveDir, 'dir')
    
    mkdir(resultSaveDir);


end

fName = ['case_', num2str(Case), '_', num2str(numel(scensToLoad)), '.mat'];
% if exist(fullfile(resultSaveDir,fName))
%     fName = [fName, '_tmp']
% 
% end
save(fullfile(resultSaveDir, [fName]),'dvhExp', 'dvhs', 'SDVH', 'SDVHdist_exp', 'SDVHdist', 'dvhDoseGrid','stdGrid','meanSDVH', '-v7.3');


%% Vis


blue = [0.35 0.7 0.9];
orange = [0.9,0.6,0];
green = [0.4660 0.6740 0.1880];
red = [1 0 0];
colors = [blue; orange; green;red];

newCT = ct;
newCT.cubeDim = dij.doseGrid.dimensions;
% newCT.cube = cellfun(@(cube) interp3(dij.ctGrid.x, dij.ctGrid.y', dij.ctGrid.z,...
%                         cube, ...
%                         dij.doseGrid.x, dij.doseGrid.y', dij.doseGrid.z), ct.cube, 'UniformOutput',false);

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
slice = 73;
profile = 153;
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
structInclude = [4,8];
structCounter = 0;
for structIdx =structInclude
    structCounter = structCounter +1;
    p(structIdx) = plot(gca(),stdGrid, SDVH(structIdx).volumePoints, 'Color',colors(structCounter,:), 'Marker', 'none', 'LineStyle','-', 'LineWidth',1.5);
    hold on;
    leg = [leg, {SDVH(structIdx).name}];
end
grid on;


legend(leg);


%% Compare results
casesToInclude = [2,4];
variables = {'dvhExp','dvhs', 'SDVH', 'SDVHdist_exp', 'SDVHdist'};
clear case_1 case_2 case_3 case_4;
clear(variables{:});
for caseIdx=casesToInclude

    load(fullfile(resultSaveDir, ['case_', num2str(caseIdx),num2str(numel(scensToLoad)), '.mat']));

    for varIdx=1:numel(variables)
        currCase.(variables{varIdx}) = eval(variables{varIdx});
       
    end
    eval(['case_', num2str(caseIdx),' = currCase;']);

    clear currCase;
end

%% compare DVHs
figure;
movegui(gcf(), 'southwest');

tiledlayout(2,1);

structInclude = [8,4];
for structIdx = structInclude
    nexttile;
    leg = [];
    caseCounter = 0;
    for caseIdx=casesToInclude
        caseCounter = caseCounter +1;
        currCase = eval(['case_', num2str(caseIdx)]);
        p(structIdx) = plot(gca(),dvhDoseGrid, currCase.dvhExp(structIdx).volumePoints, 'Color',colors(caseCounter,:), 'Marker', 'none', 'LineStyle','-', 'LineWidth',1.5);
        hold on;
    
        for i=1:size(currCase.dvhs(structIdx).volumePoints,1)
            f = plot(gca(),dvhDoseGrid, currCase.dvhs(structIdx).volumePoints(i,:), 'Color',colors(caseCounter,:), 'Marker', 'none', 'LineStyle','--', 'LineWidth',0.5);
             set(get(get(f,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
    
        leg = [leg, {['case_', num2str(caseIdx)]}];
        
    end
    grid on;
    legend(leg);
end


nexttile(1);
title('Target');
ylabel('Volume [%]');
xlabel('Dose [Gy]');

nexttile(2);
title('OAR');

ylabel('Volume [%]');
xlabel('Dose [Gy]');


%% SDVH comparison
% compare DVHs
figure;
movegui(gcf(), 'southwest');

tiledlayout(2,1);

structInclude = [8,4];
for structIdx = structInclude
    nexttile;
    leg = [];
    caseCounter = 0;
    for caseIdx=casesToInclude
        caseCounter = caseCounter +1;
        currCase = eval(['case_', num2str(caseIdx)]);
        p(structIdx) = plot(gca(),stdGrid, currCase.SDVH(structIdx).volumePoints, 'Color',colors(caseCounter,:), 'Marker', 'none', 'LineStyle','-', 'LineWidth',1.5);
        hold on;
    
        leg = [leg, {['case_', num2str(caseIdx)]}];
        
    end
    grid on;
    legend(leg);
end

nexttile(1);
title('Target');
ylabel('Volume [%]');
xlabel('Std [Gy]');

nexttile(2);
title('OAR');

ylabel('Volume [%]');
xlabel('Std [Gy]');