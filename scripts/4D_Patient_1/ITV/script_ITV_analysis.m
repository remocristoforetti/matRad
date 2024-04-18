%% load data
matRad_cfg = MatRad_Config.instance();
load('C:\r408i_data\r408i_data\CTDatasetMotion\102_HM10395_333_newCst.mat');

matRad_cfg = MatRad_Config.instance();

saveDir = fullfile(matRad_cfg.matRadRoot, 'Patient_1_analysis', 'ITV');

%% load weights
Case = 200;
load(fullfile(saveDir, ['ITV4mm_planResults_', num2str(Case), '.mat']), 'w', 'costFunctions', 'cst', 'pln');

%% analysis

RBE = 1.1;

scensToLoad = [1:30];%[1:newMultiScen.totNumScen];

[dij, newMultiScen, ~] = matRad_loadDijScenarios(ct,fullfile(saveDir, 'scen'),'all', 'none');

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
    [currDij, ~, ~] = matRad_loadDijScenarios(ct,fullfile(saveDir, 'scen'),scenIdx,[],0);
    physicalDose{scenCounter} = currDij.physicalDose{currCtIdx}*w*RBE*pln.numOfFractions;

end
    
% Get DVHs

dvhs = struct('name', [], 'volumePoints', []);
stringLenght = 0;
for scenCounter=scensToLoad
    fprintf(repmat('\b',1,stringLenght));
    stringLenght = fprintf('\t DVHs: %u/%u\n', scenCounter,numel(scensToLoad));


    %[currCtIdx, ~,~] = ind2sub(size(dij.physicalDose),scenCounter);
    tmpDVH = matRad_calcDVH(cstOnGrid,physicalDose{scenCounter},[1:10],[],dvhDoseGrid);

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
%% save

resultSaveDir = fullfile(matRad_cfg.matRadRoot, 'Patient_1_analysis', '2mm_30_patient_results');
if ~exist(resultSaveDir, 'dir')
    mkdir(resultSaveDir);


end




fName = ['case_', num2str(Case), '_', num2str(numel(scensToLoad)), '.mat'];
%fName = ['phase_case_', num2str(Case), '_', num2str(numel(scensToLoad)), '.mat'];

save(fullfile(resultSaveDir, [fName]),'dvhExp', 'dvhs', 'SDVH', 'SDVHdist_exp', 'SDVHdist', 'dvhDoseGrid','stdGrid','meanSDVH', '-v7.3');
