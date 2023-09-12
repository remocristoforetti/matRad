matRad_rc;
matRad_cfg = MatRad_Config.instance();
matRad_cfg.propOpt.defaultMaxIter = 10000;
%load 'TG119.mat'
load('C:\r408i_data\r408i_data\CTDatasetMotion\102_HM10395_333.mat');

%load('PROSTATE.mat');
%cst{3,6} = [];
ct.numOfCtScen = 10;
ct.cubeHU = ct.cubeHU(1:10);
for k=1:size(cst,1)
    cst{k,4} = cst{k,4}(1:10);
end

cst{2,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(500,0));

cst{3,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(500,0));

cst{13,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(500,0));

%% meta information for treatment plan (1) 
pln.numOfFractions  = 30;
pln.radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln.machine         = 'Generic';

% beam geometry settings
pln.propStf.bixelWidth      = 3; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.gantryAngles    = [0]; % [?] ;
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
pln.propDoseCalc.doseGrid.resolution.x = 8; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 8; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 8; % [mm]

quantityOpt  = 'physicalDose';     % options: physicalDose, effect, RBExD
%=======================================> Model check error in bioModel
modelName    = 'none';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'wcScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation

pln.multScen = matRad_multScen(ct,scenGenType);
%pln.multScen.nSamples = 5;
%pln.multScen.totNumShiftScen = 2;
%load('plnMultiScen.mat');
%% stf
stf = matRad_generateStf(ct,cst,pln);

%% cst
for voiIdx=1:size(cst,1)
    if isequal(cst{voiIdx,3}, 'TARGET')
        %if ~isempty(cst{voiIdx,6})
            for objIdx=1:size(cst{voiIdx,6},2)
                cst{voiIdx,6}{objIdx}.robustness = 'PROB';
            end
        %end
    else
        for objIdx=1:size(cst{voiIdx,6},2)

                cst{voiIdx,6}{objIdx}.robustness = 'PROB'; 

        end
    end
end

%% dij

pln.propDoseCalc.clearVoxelsForRobustness = 'none'; % none, targetOnly, oarOnly, objectivesOnly, [scenario indexes];

pln.propDoseCalc.precalcProbabilisticQuantitites = true;
dij  = matRad_calcParticleDose(ct,stf, pln,cst,0);

%% fluence opt

pln.propOpt.scen4D = 'all';
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Visualize plans

for ctIdx=1:ct.numOfCtScen
    % tmpCube = dij.physicalDose{ctIdx}*resultGUI.w;
    % robustCube{ctIdx} = reshape(tmpCube, dij.doseGrid.dimensions);
    % Nominal scenaarios should be the first 10
    robustCube{ctIdx} = resultGUI.(['physicalDose_', num2str(ctIdx)]);
end

slice = 67;
f = figure;

%movegui(f, 'southwest');
for ctIdx=1:ct.numOfCtScen  
    nCol = min([5,ct.numOfCtScen]);

    nRow = ceil(ct.numOfCtScen/nCol);
    subplot(nRow, nCol, ctIdx);
    imagesc(robustCube{ctIdx}(:,:,slice));
    matRad_plotVoiContourSlice(gca(f),cst,ct.cubeHU,ctIdx,1,3,slice);
end
