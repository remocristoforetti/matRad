matRad_rc;
matRad_cfg = MatRad_Config.instance();
load('BOXPHANTOM.mat');

SaveDirectory = 'X:\matRad_dev_varRBE_robOpt_GPU';
%% Setup
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = 0;
pln.propStf.couchAngles   = 0;
pln.propStf.bixelWidth    = 3;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;


pln.radiationMode   = 'protons';
pln.machine         = 'Generic';

pln.multScen = matRad_multScen(ct,'nomScen');

pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]
pln.propDoseCalc.calcLET = true;
stf = matRad_generateStf(ct,cst,pln);

%% Single Model Optimization
modelNames           = {'MCN'};

quantityOpt         = 'RBExD';

pln.multScen = matRad_multScen(ct,'nomScen');
pln.bioParam = matRad_bioModel(pln.radiationMode, quantityOpt, modelNames{1});


dij = matRad_calcParticleDose(ct,stf,pln,cst);

GPU =1;
resultGUI = matRad_fluenceOptimization(dij,cst,pln,GPU);

save(strcat(SaveDirectory, '\', modelNames{k}, '_BOXPHANTOM_MOD.mat'), 'resultGUI');
ParamBio = dij.TissueParameters;
ParamShift_MCN = pln.multScen;
modelParam = [{ParamBio}, {ParamShift_MCN}];
save(strcat(SaveDirectory, '\', modelNames{k}, '_BOXPHANTOM_MOD_parameters.mat'), 'modelParam');
%% Visualization
f = dir(SaveDirectory);
%name = 'doubleRob_COCW_BOXPHANTOM_MOD.mat';
name = 'MCN_BOXPHANTOM_MOD.mat';
filenames = {f(:).name};
rs_name = find(strcmp(filenames,name));
rs = load(strcat(SaveDirectory, '\', f(rs_name).name));
Plan = {rs.resultGUI.physicalDose};
nBioSamples = 1;
nSamples = 4;
x = [1:size(Plan{1},2)]*ct.resolution.x;

SingleScenarios = [];

for k = 1:nBioSamples
   SingleScenarios = [SingleScenarios, {getfield(rs.resultGUI,strcat('RBExD', '_', num2str(k)))}];
end

figure;
for k = 1:nBioSamples
   y = SingleScenarios{k}(:,80,80);
   plot(x,y,'-');

   hold on;
end

y = Plan{1}(:,80,80);
plot(x,y,'-','color', 'k');

DosePrescription = zeros(ct.cubeDim);
DosePrescription(cst{2,4}{1}) = cst{2,6}{1}.parameters{1}./pln.numOfFractions;

y = DosePrescription(:,80,80);
plot(x, y, 'o-', 'color', 'r', 'MarkerSize', 1, 'LineWidth', 0.5);

%plot(x,y,'-');


grid on;
grid minor;
%legend('MCN', 'WED', 'MKMLET', 'CAR', 'physicalDose','Dose Prescription');

%figure;
L = rs.resultGUI.LET;
y = L(:,80,80);
yyaxis right;
plot(x,y,'--');
xlabel('Depth [mm]');
ylabel('[KeV/\mu m]');
grid on;
grid minor;
%legend('MCN', 'WED', 'MKMLET', 'CAR', 'physicalDose','Dose Prescription', 'LET distribution');
%legend('MCN', 'WED', 'CAR', 'physicalDose','Dose Prescription', 'LET distribution');

legend('MCN', 'physicalDose','Dose Prescription', 'LET distribution');

%% Visualize comparison doubleRobust and ShiftRobust
f = dir(SaveDirectory);
name_double = 'doubleRob_COCW_BOXPHANTOM_MOD.mat';
%name_double = 'MCN_BOXPHANTOM_MOD_modLET.mat';
name_single = 'MCN_BOXPHANTOM_MOD.mat';
filenames = {f(:).name};
rs_name = find(strcmp(filenames,name_double));
rs_double = load(strcat(SaveDirectory, '\', f(rs_name).name));
Plan_double = {rs_double.resultGUI.physicalDose};

rs_name = find(strcmp(filenames,name_single));
rs_single = load(strcat(SaveDirectory, '\', f(rs_name).name));
Plan_single = {rs_single.resultGUI.physicalDose};

figure;
y = Plan_double{1}(:,80,80);
plot(x,y,'-');

hold on;

y = Plan_single{1}(:,80,80);

plot(x,y,'-');

grid on;
grid minor;

L_double = rs_double.resultGUI.LET;
y = L_double(:,80,80);

yyaxis right;
plot(x,y,'--', 'color', [0.9290 0.6940 0.1250]);

L_single = rs_single.resultGUI.LET;
y = L_single(:,80,80);
yyaxis right;
plot(x,y,'--', 'color', [0.4770 0.6740 0.1880]);

xlabel('Depth [mm]');
ylabel('[KeV/\mu m]');
grid on;
grid minor;
legend('physicalDose doubleRob', 'physicalDose RangeShift', 'LET doubleRob', 'LET RangeShift');