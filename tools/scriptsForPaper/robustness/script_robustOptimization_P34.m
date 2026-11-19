matRad_cfg = MatRad_Config.instance();
load('Lung_Patient_34.mat');

%% Pln
pln.numOfFractions  = 30;
pln.radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln.machine         = 'Generic';

% beam geometry settings
pln.propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.gantryAngles    = [50   90  350]; % [?] ;
pln.propStf.couchAngles     = zeros(numel(pln.propStf.gantryAngles),1); % [?] ; 
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

pln.propDoseCalc.calcLET = 0;

pln.propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln.propOpt.spatioTemp      = 0;
pln.propOpt.STscenarios     = 1;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

quantityOpt  = 'RBExD';     % options: physicalDose, effect, RBExD
%=======================================> Model check error in bioModel
modelName    = 'constRBE';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode, modelName);

pln.multScen = matRad_WorstCaseScenarios(ct);
pln.multScen.wcSigma = 2;
pln.numOfFractions = 30;
%% Stf
stf = matRad_generateStf(ct,cst,pln);
%% Dose calc
dij = matRad_calcDoseInfluence(ct,cst,stf,pln);

save(fullfile('D:\r408i\P34wc2Scenarios\dij.mat'), 'dij', 'stf', 'pln', '-v7.3');

%% cst
for i=1:size(cst,1)
    cst{i,6} = {};
    cst{i,5}.Priority = 2;
end

cst{1,5}.Priority = 1;

cst{1,6}{1} = DoseConstraints.matRad_MinMaxDose(58,65);

cst{1,6}{2} = DoseObjectives.matRad_SquaredDeviation(100,60);
cst{1,6}{2}.robustness = 'STOCH';

cst{6,6} = {DoseObjectives.matRad_MeanDose(10,0)};
cst{6,6}{1}.robustness = 'STOCH';

cst{6,6}{2} = DoseObjectives.matRad_SquaredOverdosing(100,40);
cst{6,6}{2}.robustness = 'STOCH';

cst{10,6} = {DoseObjectives.matRad_SquaredOverdosing(10,0)};
cst{10,6}{1}.robustness = 'STOCH';

%% Optimize
ct = matRad_getWorldAxes(ct);

% Distribution module
% distVisModuleProperties.quantity = 'constRBEDose';
% distVisModuleProperties.doseGrid = dij.doseGrid;
% distVisModuleProperties.ct       = ct;
% 
% distVisModuleProperties.cst      = cst;
% distVisModuleProperties.plane    = 3;
% distVisModuleProperties.slice    = 57;
% distModule           = matRad_DistributionVisModule(distVisModuleProperties);
% 
% cstOnGrid = matRad_resizeCstToGrid(cst, ct.x,ct.y,ct.z, dij.doseGrid.x, dij.doseGrid.y, dij.doseGrid.z);
% dvhModule = matRad_DVHVisualizationModule('constRBEDose', matRad_setOverlapPriorities(cstOnGrid));
% 
% % assign it to the manager
% pln.propOpt.visualizationManager.visualize = 'on';
% 
% % pln.propOpt.visualizationManager.modules = {distModule, dvhModule, matRad_ObjectiveFunctionVisModule(cst)};
% pln.propOpt.visualizationManager.modules = {matRad_ObjectiveFunctionVisModule(cst)};
matRad_cfg.defaults.propOpt.maxIter = 15000;

resultGUI = matRad_fluenceOptimization(dij,cst,pln);

savePath = fullfile(matRad_cfg.primaryUserFolder,'dataPaper', 'robustness', 'P34');
wOpt = resultGUI.wUnsequenced;
info = resultGUI.info;

if ~exist(savePath, 'dir')
    mkdir(savePath);

end
fileName = fullfile(savePath, sprintf('STOCH.mat'));

save(fileName, 'wOpt', 'pln', 'cst', 'info');
%% Vis

newCT = ct;

newCT.cubeDim    = dij.doseGrid.dimensions;
newCT.dimensions = newCT.cubeDim;

newCT.cubeHU = cellfun(@(cube) interp3(dij.ctGrid.x, dij.ctGrid.y', dij.ctGrid.z,...
                        cube, ...
                        dij.doseGrid.x, dij.doseGrid.y', dij.doseGrid.z), ct.cubeHU, 'UniformOutput',false);
newCT.cube = cellfun(@(cube) interp3(dij.ctGrid.x, dij.ctGrid.y', dij.ctGrid.z,...
                        cube, ...
                        dij.doseGrid.x, dij.doseGrid.y', dij.doseGrid.z), ct.cubeHU, 'UniformOutput',false);

newCT.x = dij.doseGrid.x;
newCT.y = dij.doseGrid.y;
newCT.z = dij.doseGrid.z;

newCT.resolution = dij.doseGrid.resolution;

newCT = rmfield(newCT, {'size', 'origin'});

cube = dij.physicalDose{1}*wOpt*dij.RBE;

cube = reshape(cube, dij.doseGrid.dimensions);
cstOnGrid = matRad_resizeCstToGrid(cst, ct.x,ct.y,ct.z, dij.doseGrid.x, dij.doseGrid.y, dij.doseGrid.z);
cstOnGrid(:,7) = cell(size(cst,1),1);

f = figure('WindowState','maximized');
% f.Position([1:2]) = f.Position([1:2]) -900;
matRad_plotSlice(newCT,'dose', cube,... 
                    'axesHandle', gca(f),...
                    'cst', cstOnGrid,...
                    'cubeIdx',1,...
                    'plane',3,...
                    'slice', 63);

%% DVH
cstOnGrid = matRad_resizeCstToGrid(cst, ct.x,ct.y,ct.z, dij.doseGrid.x, dij.doseGrid.y, dij.doseGrid.z);
dvh = matRad_calcDVH(cstOnGrid,dij.physicalDose{1}*resultGUI.wUnsequenced*dij.RBE);


figure;
for i=1:size(cst,1)
    plot(dvh(i).doseGrid, dvh(i).volumePoints, '-', 'DisplayName',strrep(cst{i,2}, '_', ' '));
    hold on;
end
grid on;

legend;
