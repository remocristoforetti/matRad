% load cst
% matRad_rc;
% load('BOXPHANTOM.mat');
% old_ct = ct;
% old_cst = cst;

%% Add OAR
% builder = matRad_PhantomBuilder(ct.cubeDim,[ct.resolution.x, ct.resolution.y, ct.resolution.z],1);
% builder.addBoxOAR('Distal_Organ',[20,10,20], 'offset', [0,20,0], 'objective', {DoseObjectives.matRad_SquaredOverdosing(400,5)});
% 
% [~, cstOAR] = builder.getctcst;
% 
% cst = [old_cst; cstOAR];
% cst{3,1} = 2;
% 

% figure;
% 
% tiledlayout(2,2);
% nexttile;
% matRad_plotVoiContourSlice(gca(),cst,ct,1,1,1,80);
% nexttile;
% matRad_plotVoiContourSlice(gca(),cst,ct,1,1,2,80);
% nexttile;
% matRad_plotVoiContourSlice(gca(),cst,ct,1,1,3,80);
% save('phantoms\BOXPHANTOM_OAR_dis.mat', 'ct', 'cst');

%% load

load('BOXPHANTOM_OAR_dis.mat');
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
pln.multScen.nSamples = 10;

stf = matRad_generateStf(ct,cst,pln);

saveDir = fullfile(matRad_cfg.matRadRoot, 'BOXPHANTOM_OAR_dist_analysis', '2mm');

[dij, dijTemplate]  = matRad_calcParticleDoseMultipleScenarios(ct,stf,pln,cst,0,saveDir,0);

%% save stf
load(fullfile(saveDir, 'stf.mat'), 'stf');