%% load ct
load('BOXPHANTOM_OAR_dis.mat');
matRad_cfg = MatRad_Config.instance();
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


%% load scenarios
saveDir = fullfile(matRad_cfg.matRadRoot, 'BOXPHANTOM_OAR_dist_analysis', '2mm');

scenMode = 'all';
[dij, newMultiScen, loadedScenarios] = matRad_loadDijScenarios(ct, saveDir, scenMode);

pln.multScen = newMultiScen;
nScens = pln.multScen.totNumScen;


%% compute Omega
cstOnGrid = matRad_resizeCstToGrid(cst,dij.ctGrid.x, dij.ctGrid.y, dij.ctGrid.z, dij.doseGrid.x, dij.doseGrid.y, dij.doseGrid.z);

dij_acc = matRad_calculateProbabilisticQuantitiesGPU(dij,cstOnGrid,pln,'all');
omega = dij_acc.physicalDoseOmega;
exp = dij_acc.physicalDoseExp;
dij.physicalDose = {[]};


save(fullfile(saveDir, 'probQuantities.mat'), 'exp', 'omega', 'dij', 'pln', 'cst', '-v7.3');