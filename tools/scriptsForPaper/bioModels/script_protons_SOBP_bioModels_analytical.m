%% Script: RBExDose SOBP in water with the analytical LET-based proton bio models
% Optimizes an SOBP on the water phantom for RBExDose with the McNamara
% model (HongPB), then recalculates the optimized weights with all the
% other LET-based proton models (and constant RBE as reference).
clear;
matRad_rc;
matRad_cfg = MatRad_Config.instance();

load(fullfile(matRad_cfg.primaryUserFolder,'dataPaper','phantomBuilder','phantom_SOBP_2.mat'),'ct','cst');

cst{1,6}{1} = DoseObjectives.matRad_SquaredDeviation(100,60);
cst{2,6}{1} = DoseObjectives.matRad_SquaredOverdosing(100,50);

pln.radiationMode   = 'protons';
pln.machine         = 'Generic';
pln.bioModel        = 'MCN';
pln.multScen        = 'nomScen';
pln.numOfFractions  = 30;

pln.propStf.gantryAngles = 270;
pln.propStf.couchAngles  = 0;
pln.propStf.bixelWidth   = 3;
pln.propStf.isoCenter    = matRad_getIsoCenter(cst,ct,0);

pln.propDoseCalc.engine = 'HongPB';
pln.propDoseCalc.doseGrid.resolution.x = 1;
pln.propDoseCalc.doseGrid.resolution.y = 1;
pln.propDoseCalc.doseGrid.resolution.z = 1;

pln.propOpt.quantityOpt = 'RBExDose';
pln.propOpt.runDAO      = false;

%% Optimize for RBExDose with McNamara
stf       = matRad_generateStf(ct,cst,pln);
dij       = matRad_calcDoseInfluence(ct,cst,stf,pln);
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
w         = resultGUI.w;

%% Recalculate the optimized weights with all models
bioModels = {'constRBE','MCN','WED','CAR','LSM'};

physicalDose = resultGUI.physicalDose;
LET          = resultGUI.LET;
RBExDose     = struct();
for i = 1:numel(bioModels)
    plnRecalc          = pln;
    plnRecalc.bioModel = bioModels{i};

    resultRecalc = matRad_calcDoseForward(ct,cst,stf,plnRecalc,w);
    RBExDose.(bioModels{i}) = resultRecalc.RBExDose;
end

outDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','bioModels');
if ~isfolder(outDir)
    mkdir(outDir);
end
save(fullfile(outDir,'protons_SOBP_bioModels_analytical_2.mat'),'ct','cst','pln','stf','w', ...
    'bioModels','physicalDose','LET','RBExDose');

%% Depth profiles along the central beam axis
isoIdx = matRad_world2cubeIndex(pln.propStf.isoCenter,ct); % [yIdx xIdx zIdx]
yIdx = isoIdx(1) + [-1,0,1];
zIdx = isoIdx(3) + [-1,0,1];
getProfile = @(cube) squeeze(mean(cube(yIdx,:,zIdx), [1,3]));

figure('Color','w');
yyaxis left;
hold on;
plot(ct.x, getProfile(physicalDose)*1.1, 'k--', 'LineWidth', 2, 'DisplayName', 'physical dose x 1.1');
for i = 1:numel(bioModels)
    plot(ct.x, getProfile(RBExDose.(bioModels{i})), 'LineWidth', 2, 'DisplayName', bioModels{i});
end
ylabel('RBExDose [Gy(RBE)]');
yyaxis right;
plot(ct.x, getProfile(LET), ':', 'LineWidth', 1.5, 'DisplayName', 'LET_d');
ylabel('LET_d [keV/\mum]');
grid on;
xlabel('x [mm]');
legend('show','Location','northwest');
