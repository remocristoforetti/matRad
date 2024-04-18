%load('4D_BOXPHANTOM_OAR_dist.mat');
% load('BOXPHANTOM_OAR_dis.mat');
% [~, cst_tmp] = matRad_addMovement(ct,cst, 5, 2, [5, 0, 0]);
% 
% ctEdge = ct;
% ctEdge.cube = ct.cube(1);
% ctEdge.cubeHU = ct.cubeHU(1);
% cstEdge = cst;
% 
% cstEdge{2,4}{1} = cst_tmp{2,4}{2};
% 
% [~, cst_mov] = matRad_addMovement(ctEdge,cstEdge,5,20,[-10,0,0]);
% 
% 
% newCt = ct;
% 
% for i=1:10
% 
%     currCube = ct.cubeHU{1};
%     currCube((cst_mov{2,4}{i})) = 350;
%     newCt.cube{i} = ct.cube{1};
%     newCt.cubeHU{i} = currCube;
% end
% newCt.numOfCtScen = 10;
% 
% 
% newCst = cst;
% 
% newCst{1,4} = repmat(cst{1,4}(1),1,10);
% newCst{2,4} = cst_mov{2,4}(1:10);
% newCst{3,4} = repmat(cst{3,4}(1),1,10);
% 
% old_ct = ct;
% ct = newCt;
% old_cst = cst;
% 
% cst = newCst;
% 
% %save('phantoms/4D_BOXPHANTOM_OAR_dist.mat', 'cst', 'ct');
% %% Vis
% figure;
% movegui(gcf(), 'southwest');
% for i=1:newCt.numOfCtScen
%     clf;
%     imagesc(ct.cubeHU{i}(:,:,80));
%     matRad_plotVoiContourSlice(gca(),cst,ct,i,[],3,80);  
%     pause(0.5);
%     i
% end
% 
% 
% 
% 
% %% save
% save('phantoms/4D_BOXPHANTOM_OAR_dist_heavy.mat', 'ct', 'cst');

%% load
load('4D_BOXPHANTOM_OAR_dist_heavy.mat');

% select only 3 scenarios
cubeIdx = [1,3,5];
cubes = [];
cubesHU = [];

for i=1:3

    cubes{i}   = ct.cube{cubeIdx(i)};
    cubesHU{i} = ct.cubeHU{cubeIdx(i)};
end

ct.cube = cubes;
ct.cubeHU = cubesHU;
ct.numOfCtScen = 3;
%% load plan parameters
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

%% Scenario setup
matRad_cfg = MatRad_Config.instance();

% stf generated on nominal scenario only
%stf = matRad_generateStf(ct,cst,pln);


% %% Run single optimization
% 
 % dij = matRad_calcParticleDose(ct,stf,pln,cst);
 % resultGUI = matRad_fluenceOptimization(dij,cst,pln);
 % w = resultGUI.w;
% save(fullfile(saveDir, 'w.mat'), 'w');
%% Compute scenarios
pln.multScen = matRad_multScen(ct,'rndScen');
pln.multScen.nSamples = 3;

stf = matRad_generateStf(ct,cst,pln, 0,[3]);


saveDir = fullfile(matRad_cfg.matRadRoot, '4D_BOXPHANTOM_OAR_dist_heavy_analysis', '2mm_3CT');


[dij, dijTemplate]  = matRad_calcParticleDoseMultipleScenarios(ct,stf,pln,cst,0,saveDir,0);

%% save stf
save(fullfile(saveDir, 'stf.mat'), 'stf');
%load(fullfile(saveDir, 'stf.mat'), 'stf');

%% check scenarios
scen1 = load(fullfile(saveDir, 'scenario_1.mat'));
scen2 = load(fullfile(saveDir, 'scenario_3.mat'));

w = ones(dij.totalNumOfBixels,1);

scen1 = scen1.dijScenario{1}*w;
scen2 = scen2.dijScenario{1}*w;

