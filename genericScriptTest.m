matRad_rc;
matRad_cfg = MatRad_Config.instance();
matRad_cfg.propOpt.defaultMaxIter = 10000;
%load 'TG119.mat'
%load('C:\r408i_data\r408i_data\CTDatasetMotion\102_HM10395_333.mat');

load('PROSTATE.mat');
%cst{3,6} = [];
%% meta information for treatment plan (1) 
pln.numOfFractions  = 5;
pln.radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln.machine         = 'Generic';

% beam geometry settings
pln.propStf.bixelWidth      = 3; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.gantryAngles    = [90 270]; % [?] ;
pln.propStf.couchAngles     = zeros(numel(pln.propStf.gantryAngles),1); % [?] ; 
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln.propDoseCalc.calcLET = 1;

pln.propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln.propOpt.spatioTemp      = 0;
pln.propOpt.STscenarios     = 2;
%pln.propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 8; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 8; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 8; % [mm]
% pln.propDoseCalc.doseGrid.resolution = ct.resolution;
quantityOpt  = 'physicalDose';     % options: physicalDose, effect, RBExD
%=======================================> Model check error in bioModel
modelName    = 'none';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'rndScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation

pln.multScen = matRad_multScen(ct,scenGenType);
pln.multScen.nSamples = 5;

load('plnMultiScen.mat');
%% stf
stf = matRad_generateStf(ct,cst,pln);

%% cst
for voiIdx=1:size(cst,1)

    for objIdx=1:size(cst{voiIdx,6},2)
        if ~isempty(cst{voiIdx,6})
            cst{voiIdx,6}{objIdx}.robustness = 'none'; 
        end
    end


end
method = 'STOCH';
cst{1,6}{1}.robustness = method;
cst{6,6}{1}.robustness = method;
cst{7,6}{1}.robustness = method;
cst{8,6}{1}.robustness = method;
%% dij
pln.propDoseCalc.clearMultiScenarioUnusedVoxels = false;

tic


dij_nominal  = matRad_calcParticleDose(ct,stf, pln,cst,0);

nominal_dij_time = toc;


tic
resultGUI_nominal = matRad_fluenceOptimization(dij_nominal,cst,pln);
nominal_opt_time = toc;



%% mod
tic
pln.propDoseCalc.clearMultiScenarioUnusedVoxels = true;

dij_reduced  = matRad_calcParticleDose(ct,stf, pln,cst,0);

reduced_dij_time = toc;


tic
resultGUI_reduced = matRad_fluenceOptimization(dij_reduced,cst,pln);
reduced_opt_time = toc;

%% opt
pln.propOpt.clearUnusedVoxels = 0;

tic
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
originalTime = toc;

pln.propOpt.clearUnusedVoxels = 1;

tic
resultGUI_reduced = matRad_fluenceOptimization(dij,cst,pln);
afterTime = toc;

%% plots
figure;
subplot(1,2,1);
imagesc(resultGUI.physicalDose(:,:,80));
title('nominal plan');
subplot(1,2,2);
reducedDistribution = matRad_calcCubes(resultGUI_reduced.w, dij,1);
imagesc(reducedDistribution.physicalDose(:,:,80));
title('reduced plan');

diff = resultGUI.physicalDose - reducedDistribution.physicalDose;
max(max(max(diff)))
w_diff = resultGUI.w - resultGUI_reduced.w;
sum(w_diff)

%% Systematic

resolutions = [8,5,3,2];
originalTime = [];
originalIterations = [];
w_var = {};
w_var_tot = [];
afterTime = [];
afterIter = [];
for resolutionIdx = 1:size(resolutions,2)
    pln.propDoseCalc.doseGrid.resolution.x = resolutions(resolutionIdx);
    pln.propDoseCalc.doseGrid.resolution.y = resolutions(resolutionIdx);
    pln.propDoseCalc.doseGrid.resolution.z = resolutions(resolutionIdx);
    stf = matRad_generateStf(ct,cst,pln);
    dij  = matRad_calcParticleDose(ct,stf, pln,cst,0);

    pln.propOpt.clearUnusedVoxels = 0;

    tic
    resultGUI = matRad_fluenceOptimization(dij,cst,pln);
    originalTime(resolutionIdx) = toc;
    originalIterations(resolutionIdx) = resultGUI.info.iter;

    pln.propOpt.clearUnusedVoxels = 1;
    
    tic
    resultGUI_reduced = matRad_fluenceOptimization(dij,cst,pln);
    afterTime(resolutionIdx) = toc;
    afterIter(resolutionIdx) = resultGUI_reduced.info.iter;

    w_var{resolutionIdx} = resultGUI.w - resultGUI_reduced.w;
    w_var_tot(resolutionIdx) = sum(w_var{resolutionIdx});
end




%% plots
figure;
plot(resolutions, originalTime, '.-');
hold on;
plot(resolutions, afterTime, '.-');
grid on;
grid minor;
legend('Original dij','Reduced dij');
xlabel('resolution [mm]', 'FontSize', 14);
ylabel('time[s]', 'FontSize', 14);
title('Prostate physicalDose', 'FontSize', 14);