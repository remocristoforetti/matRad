%% Script: FRED recalculation of the McNamara-optimized RBExDose SOBP
% Takes the plan optimized in script_protons_SOBP_bioModels_analytical.m
% and recalculates it with FRED (dose + LETd scoring) using a LET-based
% bio model. Three steps, each toggled below: write the FRED input files,
% run FRED, read the results back into matRad.

clear;
matRad_rc;
matRad_cfg = MatRad_Config.instance();

bioModel = 'MCN';

writeFiles = true;
runFred    = true;
readFred   = true;
fredUseGPU = true;


numHistoriesDirect = 1e9;

%% Load the optimized plan
anaFile = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','bioModels','protons_SOBP_bioModels_analytical_2.mat');
load(anaFile,'ct','cst','pln','stf','w','physicalDose','LET','RBExDose');

pln.bioModel = bioModel;
% the pln still selects the Hong engine used for the optimization
if isfield(pln.propDoseCalc,'engine')
    pln.propDoseCalc = rmfield(pln.propDoseCalc,'engine');
end

%% Explicit McNamara RBExDose from the loaded analytical physical dose and LETd
% alphaX/betaX cubes as matRad builds them: model defaults everywhere,
% overwritten by cst{i,5}.alphaX/betaX inside each structure
modelMCN = matRad_MCNamara();
alphaX = 0.1 * ones(ct.cubeDim);
betaX  = 0.05  * ones(ct.cubeDim);
for i = 1:size(cst,1)
    if isfield(cst{i,5},'alphaX') && isfield(cst{i,5},'betaX')
        alphaX(cst{i,4}{1}) = cst{i,5}.alphaX;
        betaX(cst{i,4}{1})  = cst{i,5}.betaX;
    end
end

abX = alphaX./betaX;

% physicalDose is per fraction, LET is the dose-averaged LETd [keV/um]
RBEmax = matRad_MCNamara.p0_MCN + matRad_MCNamara.p1_MCN * LET ./ abX;
RBEmin = matRad_MCNamara.p2_MCN + matRad_MCNamara.p3_MCN * sqrt(abX) .* LET;

effectMCN = alphaX .* RBEmax .* physicalDose + betaX .* RBEmin.^2 .* physicalDose.^2;
RBExDoseExplicit = (sqrt(alphaX.^2 + 4*betaX.*effectMCN) - alphaX) ./ (2*betaX);

fredDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','bioModels','FRED',bioModel);
if ~isfolder(fredDir)
    mkdir(fredDir);
end

%% FRED: write input files
if writeFiles
    engineFredWrite = DoseEngines.matRad_ParticleFREDEngine(pln);
    engineFredWrite.HUtable             = 'matRad_MCsquareDensity_FredMaterialConverter';
    engineFredWrite.workingDir          = fredDir;
    engineFredWrite.numHistoriesDirect  = numHistoriesDirect;
    engineFredWrite.externalCalculation = 'write';
    engineFredWrite.primaryMass = 1;
    engineFredWrite.calcDoseForward(ct,cst,stf,w);
end

%% FRED: run
if runFred
    fredRunDir = fullfile(fredDir,'MCrun');
    if ~isfile(fullfile(fredRunDir,'fred.inp'))
        matRad_cfg.dispError('No fred.inp found in %s',fredRunDir);
    end

    flags = '-f fred.inp';
    if ~fredUseGPU
        flags = ['-nogpu ' flags];
    end
    fredCall = [strtrim(DoseEngines.matRad_ParticleFREDEngine.cmdCall) ' ' flags];

    startDir = pwd;
    cd(fredRunDir);
    status = system(fredCall,'-echo');
    cd(startDir);
    if status ~= 0
        matRad_cfg.dispError('FRED exited with error code %d',status);
    end
end

%% FRED: read results
if readFred
    engineFred = DoseEngines.matRad_ParticleFREDEngine(pln);
    engineFred.externalCalculation = fredDir;
    resultFred = engineFred.calcDoseForward(ct,cst,stf,w);

    physicalDoseFred = resultFred.physicalDose;
    RBExDoseFred     = resultFred.RBExDose;
    LETFred          = resultFred.LET .* resultFred.physicalDose;
    LETDose          = LET .* physicalDose;

    save(fullfile(matRad_cfg.primaryUserFolder,'dataPaper','bioModels',sprintf('protons_SOBP_bioModels_FRED_%s.mat',bioModel)), ...
        'ct','cst','pln','stf','w','bioModel','physicalDoseFred','RBExDoseFred','LETFred');

    %% Compare with the analytical (HongPB) calculation, same model
    isoIdx = [25,1,25]; %matRad_world2cubeIndex(pln.propStf.isoCenter,ct); % [yIdx xIdx zIdx]
    yIdx = isoIdx(1) + [-5:5];% + [-1,0,1];
    zIdx = isoIdx(3) + [-5:5];% + [-1,0,1];
    getProfile = @(cube) squeeze(mean(cube(yIdx,:,zIdx), [1,3]));

    figure('Color','w');
    yyaxis left;
    hold on;
    plot(ct.x, getProfile(physicalDose),            '-',  'LineWidth', 2, 'DisplayName', 'physical dose HongPB');
    plot(ct.x, getProfile(physicalDoseFred),        '--', 'LineWidth', 2, 'DisplayName', 'physical dose FRED');
    plot(ct.x, getProfile(RBExDose.(bioModel)),     '-',  'LineWidth', 2, 'DisplayName', ['RBExDose HongPB ' bioModel]);
    plot(ct.x, getProfile(RBExDoseFred),            '--', 'LineWidth', 2, 'DisplayName', ['RBExDose FRED ' bioModel]);
    % plot(ct.x, getProfile(RBExDoseExplicit),        ':',  'LineWidth', 2, 'DisplayName', 'RBExDose explicit MCN (HongPB D, LET_d)');
    ylabel('dose [Gy] / RBExDose [Gy(RBE)]');
    yyaxis right;
    plot(ct.x, getProfile(LETDose),     ':',  'LineWidth', 1.5, 'DisplayName', 'LET_d HongPB');
    plot(ct.x, getProfile(LETFred), '-.', 'LineWidth', 1.5, 'DisplayName', 'LET_d FRED');
    ylabel('LET_d [keV/\mum]');
    grid on;
    xlabel('x [mm]');
    legend('show','Location','northwest');
end
