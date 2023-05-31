%load generic matRad environment and phantom
clear all;
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

%% Optimize the physical dose distribution for the selected plan.
% This step does not require the computation of specific biological quantities
% so the field pln.bioParam can be left unset.
dij_physicalDose = matRad_calcParticleDose(ct,stf,pln,cst,0);

pln.propOpt.quantityOpt = 'physicalDose';
resiltGUI_physicalDose = matRad_fluenceOptimization(dij_physicalDose,cst,pln);

%% Visualize the physical dose profile
Slice   = ceil(ct.cubeDim(3)/2);
Profile = ceil(ct.cubeDim(2)/2);

%define a depths vector for plotting
depths = [1:ct.cubeDim(1)]*ct.resolution.y - ct.resolution.y/2;

%plot the MCN profile for RBExD
f = figure;
subplot(1,2,1);
imagesc(resiltGUI_physicalDose.physicalDose(:,:,Slice));
matRad_plotVoiContourSlice(gca(f), cst,ct, 1, 1,3,Slice);

subplot(1,2,2);
plot(depths, resiltGUI_physicalDose.physicalDose(:,Profile,Slice), '.-');
grid on;
grid minor;
xlabel('depth [mm]', 'FontSize', 14);
ylabel(' Dose [Gy]', 'FontSize', 14);
lege = {'physicalDose'};
legend(lege, 'FontSize', 14);

%% Apply model for biological optimization
% In this case, a specific biological model needs to be selected. For
% example, choose the McNamara LET based model for protons 

modelName = 'MCN';
pln.bioParam = matRad_bioModel(pln.radiationMode,modelName);

%% Compute the dij for the MCN model
% With respect to the previously computed dij for physical dose only, this
% struct will also contain the model-specific parameters used for RBE and 'effect' computation 
dij_MCN = matRad_calcParticleDose(ct,stf,pln,cst,0);

%% Select optimization quantity
%We can now choose to optimize the RBE weighted dose distribution (RBExD)
%or the 'effect' distribution
pln.propOpt.quantityOpt = 'RBExD';

[resultGUI_MCN] = matRad_fluenceOptimization(dij_MCN,cst,pln);

%% Visualize the RBExD distribution and the corresponding physical dose
f = figure;
subplot(1,2,1);
imagesc(resultGUI_MCN.RBExD(:,:,Slice));
matRad_plotVoiContourSlice(gca(f), cst,ct, 1, 1,3,Slice);

subplot(1,2,2);
plot(depths, resultGUI_MCN.RBExD(:,Profile,Slice), '.-');
hold on;
plot(depths, resultGUI_MCN.physicalDose(:,Profile,Slice), '.-', 'color', 'k');

grid on;
grid minor;
xlabel('depth [mm]', 'FontSize', 14);
ylabel(' RBE weighted Dose [Gy]', 'FontSize', 14);
lege = {'RBExD', 'Physical Dose'};
legend(lege, 'FontSize', 14);

%% Other biological models can be applied to the already computed MCN plan
% This step allows for the comparison of different RBE models applied to
% the same dose distribution

modelNames = {'WED', 'CAR', 'constRBE'};

RBExD_LETmodel = [];

%loop over all the selected models
for k=1:size(modelNames,2)
    %select the model
    pln.bioParam = matRad_bioModel(pln.radiationMode,modelNames{k});
    
    %compute the dij
    dij = matRad_calcParticleDose(ct,stf,pln,cst,0);
    
    %compute the dose distribution applying the weights optimized for the
    %MCN plan
    resultGUI_LETmodel = matRad_calcCubes(resultGUI_MCN.w,dij,1);
    RBExD_LETmodel = [RBExD_LETmodel,{resultGUI_LETmodel.RBExD}];
end
%% Plot the comparison between different models

f = figure;
subplot(1,2,1);
imagesc(resultGUI_MCN.RBExD(:,:,Slice));
matRad_plotVoiContourSlice(gca(f), cst,ct, 1, 1,3,Slice);

subplot(1,2,2);
plot(depths, resultGUI_MCN.RBExD(:,Profile,Slice), '.-');
lege = {'RBExD MCN (Optimized)'};
hold on;
for k=1:size(RBExD_LETmodel,2)
    plot(depths, RBExD_LETmodel{k}(:,Profile,Slice), '.-');
    lege = [lege, {['RBExD ', modelNames{k}]}];
end

plot(depths, resultGUI_MCN.physicalDose(:,Profile,Slice), '.-', 'color', 'k');
lege = [lege, {'physicalDose'}];
grid on;
grid minor;
xlabel('depth [mm]', 'FontSize', 14);
ylabel(' RBE weighted Dose [Gy]', 'FontSize', 14);
legend(lege, 'FontSize', 14);

