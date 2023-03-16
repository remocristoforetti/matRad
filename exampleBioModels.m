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

pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]

%% Select radiation mode
pln.radiationMode   = 'protons'; % 'carbon', 'helium'

%% Select the base data to be used.
%Only machine that contains info on spectra for tabulated RBE models is
%newGeneric, for the time being. For protons, this is obtained starting
%from the generic_TOPAS_mcSquare.mat basedata. This base data can be easily
%reproduced in topas (for what concerns the range) but misses the LET
%quantity, which is needed for the LET based biolModels. So ProtonLET is
%also scored in topas when building the fluence spectra. This does not exactly 
%match the LET that is stored in generic machine for protons, need to be checked.

%for carbons, newGeneric is instead simulated starting from generic machine
%(which has LET), and an arbitrary 0.5 energy spread is added for all
%energies and seems to reproduce data quite well.
pln.machine         = 'newGeneric';

%% Generate the stf file
stf = matRad_generateStf(ct,cst,pln);

%% Override the stf information in order to select a single bixel.
%This has to be changed
%Just skip this passage if info for all the energies/full SOBP is required.
stf = OverrideStf(stf,1); %Override and keep 1 energy

%Select which energy, thi
eIdx = 59;
stf.ray.energy = machine.data(eIdx).energy;

%% Minor plan setup info
pln.propDoseCalc.calcLET = 1;

load([pln.radiationMode, '_',pln.machine]);

quantityOpt = 'RBExD';
%% Compute example of LET based bioModel
%Set the biological parameters;
modelName = 'MCN';
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName,pln.machine);

dij_LETBased_model = [];
resultGUI_LETBased_model = [];
%Compute the dij and the RBExD cubes
dij_LETBased_model = matRad_calcParticleDose(ct,stf,pln,cst,0);
resultGUI_LETBased_model = matRad_calcCubes(1,dij_LETBased_model,1);

%% Repeat calculation for other models
modelNames = {'WED', 'CAR', 'MKMLET'};

for modelIdx=1:size(modelNames,2)
    pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelNames{modelIdx},pln.machine);
    dij_currModel = matRad_calcParticleDose(ct,stf,pln,cst,0);
    currResultGUI = matRad_calcCubes(1,dij_currModel,1);

    dij_LETBased_model = [dij_LETBased_model,dij_currModel];
    resultGUI_LETBased_model = [resultGUI_LETBased_model, currResultGUI];
end

%% Plot central profile
Slice   = ct.cubeDim(3)/2;
Profile = ct.cubeDim(2)/2;

figure;
depths = [1:ct.cubeDim(1)]*ct.resolution.y - ct.resolution.y/2;

for modelIdx = 1:length(resultGUI_LETBased_model)
    profileRBExD = resultGUI_LETBased_model(modelIdx).RBExD(:,Profile,Slice);
    plot(depths, profileRBExD, '.-');
    hold on;
end

grid on;
grid minor;
xlabel('depth [mm]');

ylabel('Dose [Gy]');
legend([modelName, modelNames]);

%% Optimize a SOBP, given one specific model
stf = matRad_generateStf(ct,cst,pln);
modelName = 'MCN';
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName,pln.machine);

dij_LETBased_model = [];
resultGUI_LETBased_model = [];
%Compute the dij and the RBExD cubes
dij_LETBased_model = matRad_calcParticleDose(ct,stf,pln,cst,0);
resultGUI_LETBased_model = matRad_fluenceOptimization(dij_LETBased_model,cst,pln);

%% Recompute the same dose distribution with other RBE models

modelNames = {'WED', 'CAR', 'MKMLET'};


calcResultGUI = [];

for modelIdx=1:size(modelNames,2)
    pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelNames{modelIdx},pln.machine);

    dij_currModel = matRad_calcParticleDose(ct,stf,pln,cst,0);
    currResultGUI = matRad_calcCubes(resultGUI_LETBased_model(1).w,dij_currModel,1);
    dij_LETBased_model = [dij_LETBased_model,dij_currModel];
    calcResultGUI = [calcResultGUI,currResultGUI];

end
%% Plot
Slice   = ct.cubeDim(3)/2;
Profile = ct.cubeDim(2)/2;

figure;
depths = [1:ct.cubeDim(1)]*ct.resolution.y - ct.resolution.y/2;
profileRBExD = resultGUI_LETBased_model.RBExD(:,Profile,Slice);
plot(depths, profileRBExD, '.-');
hold on;

for modelIdx = 1:length(calcResultGUI)
    profileRBExD = calcResultGUI(modelIdx).RBExD(:,Profile,Slice);
    plot(depths, profileRBExD, '.-');
end

grid on;
grid minor;
xlabel('depth [mm]');

ylabel('Dose [Gy]');
legend([modelName, modelNames]);