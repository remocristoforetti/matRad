matRad_rc;
matRad_cfg = MatRad_Config.instance();
matRad_cfg.propOpt.defaultMaxIter = 10000;

load('C:\r408i_data\r408i_data\CTDatasetMotion\102_HM10395_333_newCst.mat');

ct.numOfCtScen = 10;
ct.cubeHU = ct.cubeHU(1:ct.numOfCtScen);
for k=1:size(cst,1)
    cst{k,4} = cst{k,4}(1:ct.numOfCtScen);
end

% Lungs
cst{2,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(500,20));
cst{3,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(500,20));

% Heart

cst{4,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(500,20));
cst{13,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,20));
%% Plan setup
pln.numOfFractions  = 30;
pln.radiationMode   = 'photons';           % either photons / protons / helium / carbon
pln.machine         = 'Generic'; %'HITfixedBL';

% beam geometry settings
pln.propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
anglesLinear = floor(linspace(0,360,8));
pln.propStf.gantryAngles    = anglesLinear(1:end-1); % [?] ;
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
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

quantityOpt  = 'physicalDose';     % options: physicalDose, effect, RBExD
%=======================================> Model check error in bioModel
modelName    = 'none';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions

%pln.multScen = matRad_multScen(ct, 'nomScen');
 scenGenType  = 'rndScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          
% 
% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

pln.multScen = matRad_multScen(ct,scenGenType);
% pln.multScen.includeNominalScenario = true;
pln.multScen.nSamples = 10; 


%% stf
stf = matRad_generateStf(ct,cst,pln);



%% dose calc
pln.propDoseCalc.precalcProbabilisticQuantitites = false;
pln.propDoseCalc.accumulateQuantities = true;
pln.propDoseCalc.probabilisticQuantitiesMode = 'all';
pln.propDoseCalc.useGPUtoAccumulateQuantitites = true;
tic
dij_accumulate  = matRad_calcPhotonDoseAccumulate(ct,stf, pln,cst,0);
acc_time = toc;
%% Usual
pln.propDoseCalc.precalcProbabilisticQuantitites = true;
pln.propDoseCalc.accumulateQuantities = false;
tic
dij = matRad_calcPhotonDose(ct,stf, pln,cst,0);
time_omegaNominal_scenariosOnly = toc;


cst_n = cst;
for i=1:size(cst,1)
    for j=1:size(cst{i,6})
        cst_n{i,6}{j} = matRad_DoseOptimizationFunction.createInstanceFromStruct(cst{i,6}{j});
    end
end

cstOnGrid = matRad_resizeCstToGrid(cst_n,dij.ctGrid.x,dij.ctGrid.y,dij.ctGrid.z,dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);

dij_omegaNominal = matRad_calculateProbabilisticQuantitiesGPU(dij,cstOnGrid,pln,'all');
%% diff

ctIdx=1;
structIdx = 8;
wOnes = ones(dij.totalNumOfBixels,1);
Exp_nominal = dij_omegaNominal.physicalDoseExp{ctIdx}*wOnes;

exp_acc = dij_accumulate.physicalDoseExp{ctIdx}*wOnes;
diff_exp = zeros(size(exp_acc));

voxIdx = unique([find(Exp_nominal>0); find(exp_acc >0)]);
diff_exp(voxIdx) = 100*(Exp_nominal(voxIdx) - exp_acc(voxIdx))./(Exp_nominal(voxIdx));

Exp_nominal_mat = reshape(Exp_nominal, dij_accumulate.doseGrid.dimensions);
exp_acc_mat = reshape(exp_acc, dij_accumulate.doseGrid.dimensions);
diff_exp_mat = reshape(diff_exp, dij_accumulate.doseGrid.dimensions);

diff_Omega = zeros(size(dij_accumulate.physicalDoseOmega{structIdx,ctIdx}));
voxIdx_omega = unique([find(dij_accumulate.physicalDoseOmega{structIdx,ctIdx}>0); find(dij_omegaNominal.physicalDoseOmega{structIdx,ctIdx} >0)]);
diff_Omega(voxIdx_omega) = (dij_accumulate.physicalDoseOmega{structIdx,ctIdx}(voxIdx_omega) - dij_omegaNominal.physicalDoseOmega{structIdx,ctIdx}(voxIdx_omega))./dij_omegaNominal.physicalDoseOmega{structIdx,ctIdx}(voxIdx_omega);

 % Vis

slice = 25;
figure;
subplot(1,3,1);
imagesc(Exp_nominal_mat(:,:,slice));
colorbar;
subplot(1,3,2);
imagesc(exp_acc_mat(:,:,slice));
colorbar;
subplot(1,3,3);
imagesc(diff_exp_mat(:,:,slice));
colorbar;
%movegui(gcf,'southwest');
colorbar;


%% Omega vis
figure;
subplot(1,3,1);
imagesc(dij_omegaNominal.physicalDoseOmega{structIdx, ctIdx});
colorbar;
subplot(1,3,2);
imagesc(dij_accumulate.physicalDoseOmega{structIdx, ctIdx});
colorbar;
subplot(1,3,3);
imagesc(diff_Omega);
colorbar;

%movegui(gcf, 'southwest');
figure;
histogram(diff_Omega);
set(gca(), 'yscale', 'log');

fprintf('min diff_omega = %d, max diff_omega = %d\n', min(diff_Omega, [], 'all'), max(diff_Omega, [], 'all'));
%% Fluence Opt
for i=1:size(cst,1)
    for j=1:size(cst{i,6},2)
        cst{i,6}{j}.robustness = 'PROB';
    end
end

includedStructs = [2,3,4,8,13];

for v=includedStructs
    cst{v,6}{2} = OmegaObjectives.matRad_TotalVariance();
    cst{v,6}{2}.penalty = cst{v,6}{1}.penalty;

end
pln.propOpt.scen4D = [];
resultGUI_acc = matRad_fluenceOptimization(dij_accumulate,cst,pln);

%% Fluence nominal

pln.propOpt.scen4D = [];
resultGUI_omegaNominal = matRad_fluenceOptimization(dij_omegaNominal,cst,pln);

fakeDij = dij_omegaNominal;
fakeDij.physicalDose{1} = dij_omegaNominal.physicalDoseExp{1};

nominalExpDist = matRad_calcCubes(resultGUI_omegaNominal.w,fakeDij,1);

%% Plans
slice = 71;
figure;
subplot(1,3,1);
matRad_plotSliceWrapper(gca(), ct,cst,1,nominalExpDist.physicalDose,3,slice);
subplot(1,3,2);
matRad_plotSliceWrapper(gca(), ct, cst,1,resultGUI_acc.physicalDoseExp,3,slice);
subplot(1,3,3);
resolutionC = [ct.resolution.x, ct.resolution.y, ct.resolution.z];
gammaC = matRad_gammaIndex(nominalExpDist.physicalDose,resultGUI_acc.physicalDoseExp,resolutionC,[3 3]);
map = gammaIndex;

matRad_plotSliceWrapper(gca(), ct,cst,1,gammaC,3,slice,[],[],[],map);
clim([0,2]);
colormap(map);
colorbar;
movegui(gcf, 'southwest');
%% Weights
figure;
plot([1:dij_accumulate.totalNumOfBixels], resultGUI_acc.w, '.-');
hold on;
plot([1:dij_omegaNominal.totalNumOfBixels], resultGUI_omegaNominal.w, '.-');

movegui(gcf, 'southwest');

