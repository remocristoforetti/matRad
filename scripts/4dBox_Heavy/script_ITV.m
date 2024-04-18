matRad_cfg = MatRad_Config.instance();
load('4D_BOXPHANTOM_OAR_dist_heavy.mat')
old_cst = cst;
old_ct = ct;


%% Get voxels
allVoxels = arrayfun(@(scenStruct) scenStruct{1}', cst{2,4}, 'UniformOutput',false);
structVoxels = unique([allVoxels{:}])';
HUTarget = ct.cubeHU{1}(allVoxels{1}(100));
%% Define ITV
itvSt{1} = [3];
itvSt{2} = 'ITV';
itvSt{3} = 'TARGET';
itvSt{4} = cell(1,ct.numOfCtScen);
for scenIdx=1:ct.numOfCtScen
    itvSt{4}{scenIdx} = structVoxels;
end
itvSt{5} = cst{2,5};
itvSt{5}.visibleColor = [0.9, 0.1, 0];
itvSt{6} = [];

cst(4,:) = itvSt;

% Density override
for scenIdx=1:ct.numOfCtScen
    ct.cubeHU{scenIdx}(structVoxels) = HUTarget;
end

%% ITV margin
marginVolume = 2; %mm

itvMargin.x = marginVolume; % mm
itvMargin.y = marginVolume; % mm
itvMargin.z = marginVolume; % mm

itvMask = zeros(ct.cubeDim);
itvMask(cst{4,4}{1}) = 1;
expItvMask = matRad_addMargin(itvMask,cst,ct.resolution,itvMargin);


itvMarginSt{1} = 4;
itvMarginSt{2} = 'ITV_2mm_margin';
itvMarginSt{3} = 'TARGET';
itvMarginSt{4} = repmat({find(expItvMask)}, 1, ct.numOfCtScen);
itvMarginSt{5} = cst{4,5};
itvMarginSt{5}.visibleColor = [0.7, 0.2, 0.1];
itvMarginSt{6} = cst{2,6};

cst(5,:) = itvMarginSt;

%% Plan setup
pln.numOfFractions  = 30;
pln.radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln.machine         = 'Generic';

% beam geometry settings
pln.propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.gantryAngles    = 0; % [?] ;
pln.propStf.couchAngles     = zeros(numel(pln.propStf.gantryAngles),1); % [?] ; 
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln.propDoseCalc.calcLET = 0;

pln.propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln.propOpt.spatioTemp      = 0;
pln.propOpt.STscenarios     = 1;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 2; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 2; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 2; % [mm]

quantityOpt  = 'RBExD';     % options: physicalDose, effect, RBExD
%=======================================> Model check error in bioModel
modelName    = 'constRBE';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

pln.multScen = matRad_multScen(ct,'nomScen');

%% stf

stf = matRad_generateStf(ct,cst,pln);

%% dose calc
dij = matRad_calcParticleDose(ct,stf,pln,cst);

%% Load dij


load(fullfile(saveDir, 'ITV_2mm.mat'), 'dij');

%% Optimization parameters

script_OptimizationParameters_PROB_4Dbox;

cst{1,6}{1}.robustness = 'none';
cst{2,6}{1}.robustness = 'none';
cst{3,6}{1}.robustness = 'none';

cst{5,6} = cst{2,6};

cst{2,6} = {};
cst{4,6} = {};
Case =2;
switch Case
    case 1
        cst{1,6} = {};
        cst{3,6} = {};
       
    case 2
        
end


%% Optimize
pln.propOpt.scen4D = 1;
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
matRad_cfg = MatRad_Config.instance();
%% save
w = resultGUI.w;
prob_info = resultGUI.info;
costFunctions = resultGUI.costFunctions;


saveDir = fullfile(matRad_cfg.matRadRoot, '4D_BOXPHANTOM_OAR_dist_heavy_analysis', 'ITV');

if ~exist(saveDir, 'dir')

    mkdir(saveDir);
end

save(fullfile(saveDir, 'ITV_2mm_parameters.mat'), 'dij', 'stf', 'cst','-v7.3');
save(fullfile(saveDir, ['plan_results_', num2str(Case),'.mat']),  'w', 'prob_info', 'costFunctions');


%% analysis
%recover the scenario shifts and probabilities
oldScenDir = fullfile(matRad_cfg.matRadRoot, '4D_BOXPHANTOM_OAR_dist_heavy_analysis','2mm');

[~,newMultiScen, ~] = matRad_loadDijScenarios(old_ct,oldScenDir,'all', 'none');

pln.multScen = newMultiScen;
%dij_analysis = matRad_calcParticleDoseMultipleScenarios(old_ct,stf,pln,cst,0,fullfile(saveDir,'scenarios'),0);

dij_analysis = matRad_loadDijScenarios(old_ct, oldScenDir, 'all');
%% laoding

Case =2;
load(fullfile(saveDir, ['plan_results_', num2str(Case), '.mat']), 'w', 'costFunctions');
%load(fullfile(saveDir, ['planResults_phase_', num2str(Case), '.mat']), 'w', 'costFunctions');


RBE = 1.1;
scensToLoad = [1:30];%[1:newMultiScen.totNumScen];

[dij, newMultiScen, ~] = matRad_loadDijScenarios(ct,fullfile(saveDir,'scenarios'),'all', 'none');

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
    [currDij, ~, ~] = matRad_loadDijScenarios(ct,fullfile(saveDir,'scenarios'),scenIdx,[],0);
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
structVoxels = [];
for structIdx=1:size(cstOnGrid,1)
    allVoxels = arrayfun(@(scenStruct) scenStruct{1}', cstOnGrid{structIdx,4}, 'UniformOutput',false);
    structVoxels{structIdx} = unique([allVoxels{:}])';

    %this is equal to vTot/#voxels -> vTot/#Voxels = mean((SDVHdist(voxels)/pln.numOfFractions)^2)
    %meanSDVH{structIdx} = mean((SDVHdist(structVoxels{structIdx})/pln.numOfFractions).^2);

    meanSDVH{structIdx} = mean(SDVHdist(structVoxels{structIdx}));
end

%% save

resultSaveDir = fullfile('4D_BOXPHANTOM_OAR_dist_heavy_analysis', '2mm_results');
if ~exist(resultSaveDir, 'dir')
    mkdir(resultSaveDir);


end

fName = ['ITV_case_', num2str(Case), '_', num2str(numel(scensToLoad)), '.mat'];

save(fullfile(resultSaveDir, [fName]),'dvhExp', 'dvhs', 'SDVH', 'SDVHdist_exp', 'SDVHdist', 'dvhDoseGrid','stdGrid','meanSDVH', '-v7.3');
