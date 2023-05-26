%load generic matRad environment and phantom
matRad_rc;
matRad_cfg = MatRad_Config.instance();
load('BOXPHANTOM.mat');
%% Setup

%Define generic plan info
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = 0;
pln.propStf.couchAngles   = 0;
pln.propStf.bixelWidth    = 3;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

pln.multScen = matRad_multScen(ct,'nomScen');

pln.propDoseCalc.doseGrid.resolution.x = 8; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 8; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 8; % [mm]

%% Select radiation mode and machine
pln.radiationMode   = 'protons'; % 'carbon', 'helium'

pln.machine         = 'Generic';

%% Generate the stf file
stf = matRad_generateStf(ct,cst,pln);

%% Select generic plan info
pln.propDoseCalc.calcLET = 1;
quantityOpt = 'RBExD';

%% Compute example of LET based bioModel
%Set the biological parameters;
modelName = 'MCN';
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName,pln.machine);

%% Compute the dij for the MCN model
dij_MCN = matRad_calcParticleDose(ct,stf,pln,cst,0);

%% Optimize the SOBP
resultGUI_MCN = matRad_fluenceOptimization(dij_MCN,cst,pln);

%% Repeat calculation for other models
%Select the models
modelNames = {'WED', 'CAR', 'LSM'};

dij_LETBased_model = [];
resultGUI_LETBased_model = [];

for modelIdx=1:size(modelNames,2)
    %Load the specific bioModel
    pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelNames{modelIdx},pln.machine);
    
    %Compute the corresponding dij
    dij_currModel = matRad_calcParticleDose(ct,stf,pln,cst,0);
    
    %Apply the weights optimized for the MCN model
    currResultGUI = matRad_calcCubes(resultGUI_MCN.w, dij_currModel,1);

    %Store the results
    dij_LETBased_model = [dij_LETBased_model,dij_currModel];
    resultGUI_LETBased_model = [resultGUI_LETBased_model, currResultGUI];

end

%% Plot central profile for MCN model
%Select a central profile
Slice   = ceil(ct.cubeDim(3)/2);
Profile = ceil(ct.cubeDim(2)/2);

%define a depths vector for plotting
depths = [1:ct.cubeDim(1)]*ct.resolution.y - ct.resolution.y/2;

%plot the MCN profile for RBExD
figure;
plot(depths, resultGUI_MCN.RBExD(:,Profile,Slice), '.-');
grid on;
grid minor;
xlabel('depth [mm]', 'FontSize', 14);
ylabel('RBE weighted Dose [Gy]', 'FontSize', 14);
lege = {'MCN'};

legend(lege, 'FontSize', 14);

%% Plot also the other models
hold on;
for modelIdx = 1:length(resultGUI_LETBased_model)
    profileRBExD = resultGUI_LETBased_model(modelIdx).RBExD(:,Profile,Slice);
    plot(depths, profileRBExD, '.-');
    lege(modelIdx+1) = modelNames(modelIdx);
end
plot(depths, resultGUI_MCN.physicalDose(:,Profile, Slice), '--', 'color', 'k');
lege = [lege, 'physical Dose'];

%plot(depths, resultGUI_MCN.LET(:,Profile,Slice), '--');
%lege = [lege, 'physical Dose', 'LET'];
legend(lege, 'FontSize', 14);

%% Repeat same procedure for effect
quantityOpt = 'RBExD';

%Set the biological parameters;
modelName = 'MCN';
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName,pln.machine);

dij_MCN_effect = matRad_calcParticleDose(ct,stf,pln,cst,0);

resultGUI_MCN_effect = matRad_fluenceOptimization(dij_MCN_effect,cst,pln);

%Select the models
modelNames = {'WED', 'CAR'};

dij_LETBased_model_effect = [];
resultGUI_LETBased_model_effect = [];

for modelIdx=1:size(modelNames,2)
    %Load the specific bioModel
    pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelNames{modelIdx},pln.machine);
    
    %Compute the corresponding dij
    dij_currModel = matRad_calcParticleDose(ct,stf,pln,cst,0);
    
    %Apply the weights optimized for the MCN model
    currResultGUI = matRad_calcCubes(resultGUI_MCN.w, dij_currModel,1);

    %Store the results
    dij_LETBased_model_effect = [dij_LETBased_model,dij_currModel];
    resultGUI_LETBased_model_effect = [resultGUI_LETBased_model_effect, currResultGUI];

end

%% Plot the effect optimization
Slice   = ceil(ct.cubeDim(3)/2);
Profile = ceil(ct.cubeDim(2)/2);

%define a depths vector for plotting
depths = [1:ct.cubeDim(1)]*ct.resolution.y - ct.resolution.y/2;

%plot the MCN profile for RBExD
figure;
plot(depths, resultGUI_MCN_effect.effect(:,Profile,Slice), '.-');
grid on;
grid minor;
xlabel('depth [mm]', 'FontSize', 14);
ylabel('RBE weighted Dose [Gy]', 'FontSize', 14);
lege = {'MCN'};

legend(lege, 'FontSize', 14);

% Plot also the other models
hold on;
for modelIdx = 1:length(resultGUI_LETBased_model_effect)
    profileRBExD = resultGUI_LETBased_model_effect(modelIdx).effect(:,Profile,Slice);
    plot(depths, profileRBExD, '.-');
    lege(modelIdx+1) = modelNames(modelIdx);
end

legend(lege, 'FontSize', 14);

%% Other bioModels
%Other biological models are available in matRad. For example, tissue-specific tabulated
%values for alpha and beta can be provided and coupled with the spectral
%fluence information stored in base Data. This allows the estimation of
%variable alpha/beta values for the mixed field radiation.
%The "RBE tables" (namely energy dependent alpha/beta for every ion), can
%be independently computed according to different models,
%approximations and parameters. Currently available are MKM and LEM (I, II, III) models.

