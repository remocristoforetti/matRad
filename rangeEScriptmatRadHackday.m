%%% Test E-R fit
%% Proton fit
clear all;
%Fit arbitrary ranges
 meanEnergy = @(x) 5.762374661332111e-20 * x^9 - 9.645413625310569e-17 * x^8 + 7.073049219034644e-14 * x^7 ...
                        - 2.992344292008054e-11 * x^6 + 8.104111934547256e-09 * x^5 - 1.477860913846939e-06 * x^4 ...
                        + 1.873625800704108e-04 * x^3 - 1.739424343114980e-02 * x^2 + 1.743224692623838e+00 * x ...
                        + 1.827112816899668e+01;
ranges = linspace(0,300,1000);
fittedEnergies = arrayfun(@(r80) meanEnergy(r80),ranges);

figure;
subplot(2,1,1);
plot(fittedEnergies, ranges, '.-');
grid on;
xlabel('E [MeV]');
ylabel('Range [mm]');

Bortfeld_alpha = 0.022;
Bortfeld_p = 1.77;

BortfeldRanges = Bortfeld_alpha.*(fittedEnergies.^Bortfeld_p);
hold on;
plot(fittedEnergies,BortfeldRanges, '.-');
legend('meanEnergy fit', 'Bortfeld fit');
title('Proton meanEnergy fit vs Bortfeld fit');
subplot(2,1,2);
plot(fittedEnergies,BortfeldRanges-ranges, '.-');
grid on;
xlabel('E [MeV]');
ylabel('Bortfeld - fitted [mm]');

%% load base data
load('protons_Generic.mat');

baseDataEnergies = [machine.data(:).energy];
baseDataRanges   = [machine.data(:).range];

%% PROTONS: BaseData range vs baseData Z r80
baseDataRangesZ = [];
%%%% MCemittance r80 sampling
for k=1:size(baseDataEnergies,2)
   newDepths = linspace(0,machine.data(k).depths(end),numel(machine.data(k).depths) * 100);
   newDose   = interp1(machine.data(k).depths, machine.data(k).Z, newDepths, 'spline');

   [maxV, maxI] = max(newDose);
   [~, r80ind] = min(abs(newDose(maxI:end) - 0.8 * maxV));
   r80ind = r80ind - 1;
   baseDataRangesZ(k) = interp1(newDose(maxI + r80ind - 1:maxI + r80ind + 1), ...
             newDepths(maxI + r80ind - 1:maxI + r80ind + 1), 0.8 * maxV);
end
figure;
subplot(2,1,1);
plot(baseDataEnergies,baseDataRanges, 'o');
hold on;
plot(baseDataEnergies,baseDataRangesZ, '.-');
grid on;
xlabel('Energy [MeV]');
ylabel('Range [mm]');
legend('baseData range', 'interpolated r80');
title('proton baseData vs interp r80');
subplot(2,1,2);
plot(baseDataEnergies,baseDataRanges-baseDataRangesZ, '.-');
grid on;
xlabel('Energy [MeV]');
ylabel('baseData - r80 [mm]');

%% PROTONS: baseData ranges vs fitted Ranges and enegies
figure;
subplot(3,1,1);
plot(baseDataEnergies,baseDataRanges, '.-');
fittedEnergies = arrayfun(@(r80) meanEnergy(r80),baseDataRanges);
fittedRangesInterp = interp1(fittedEnergies,baseDataRanges,baseDataEnergies, 'spline');
hold on;
plot(baseDataEnergies,fittedRangesInterp, '.-');
grid on;
legend('Base data', 'meanEnergy fit');
xlabel('Energy [MeV]');
ylabel('Range [mm]');
title('Protons BaseData vs Fit');

subplot(3,1,2);
plot(baseDataEnergies,baseDataRanges-fittedRangesInterp, '.-');
xlabel('Energy [MeV]');
ylabel('baseData - fitted [mm]');
grid on;

subplot(3,1,3);
plot(baseDataRanges,baseDataEnergies-fittedEnergies, '.-');
grid on;
xlabel('Range [mm]');
ylabel('baseData - fitted [MeV]');

%% PROTONS: Bortfeld vs base data
figure;
BortfeldBaseDataRanges = Bortfeld_alpha.*(baseDataEnergies.^Bortfeld_p);
subplot(2,1,1);
plot(baseDataEnergies, baseDataRanges, '.-');
hold on;
plot(baseDataEnergies, BortfeldBaseDataRanges, '.-');
grid on;
xlabel('Energy [MeV]');
ylabel('Range [mm]');
legend('baseData', 'Bortfeld');
title('Protons Bortfeld vs baseData');
subplot(2,1,2);
plot(baseDataEnergies,BortfeldBaseDataRanges-baseDataRanges, '.-');
grid on;
xlabel('Energy [MeV]');
ylabel('Bortfeld - base data [mm]');

%% PROTONS: topas check
%eIdx = [7,28,50];
eIdx = [1:10:size(machine.data,2)];
eIdx = eIdx(1:end-2);
energies = baseDataEnergies(eIdx);
wDir = 'D:\matRad_gitHubRemo_MKM_TOPAS\topas\MCrun\BaseData_protons\Results\generic_MCenergy_TOPAS_vacuum';

name = 'score_BaseData_field1_run1_physicalDose_MCenergy_vacuum_Run_';

%figure;
x = [1:1500]*0.2 - 0.1;
integration = [45:55];
%Load vaccum profiles
for k=0:size(eIdx,2)-1
   fileName = strcat(wDir, filesep,name, string(compose('%04d', k)));
   data = squeeze(double(dicomread(fileName)));
   data = permute(data, [2,3,1]);

   profileData = squeeze(sum(data(integration,integration,:),[1,2]));
   
   [maxV, maxI] = max(profileData);
   [~, r80ind] = min(abs(profileData(maxI:end) - 0.8 * maxV));
   r80ind = r80ind - 1;
   topasR80(k+1) = interp1(profileData(maxI + r80ind - 1:maxI + r80ind + 1), ...
             x(maxI + r80ind - 1:maxI + r80ind + 1), 0.8 * maxV);
   
   PDD(:,k+1) = profileData;
   %AUC_PDD(k) = sum(profileData).*(x(2) -x(1));
%    plot(x, profileData./max(profileData), '.-', 'color','r');
%    hold on;

end

%Plot Vacuum profiles
figure;
lege = [];
for k=3:size(eIdx,2)-3
   plot(x, PDD(:,k)./max(PDD(:,k)), '.-');
   hold on;
   plot(machine.data(eIdx(k)).depths,machine.data(eIdx(k)).Z./max(machine.data(eIdx(k)).Z), '--', 'color','k');
   xline(baseDataRanges(eIdx(k)));
   lege = [lege, {['Topas : ', num2str(energies(k))]}, {['baseData : ', num2str(energies(k))]}];
end
xlabel('depth [mm]');
ylabel('dose');
legend(lege);
grid on;
title('Vacuum');

%Load air profiles
wDir = 'D:\matRad_gitHubRemo_MKM_TOPAS\topas\MCrun\BaseData_protons\Results\generic_MCenergy_TOPAS_air';
name = 'score_BaseData_field1_run1_physicalDose_MCenergy_vacuum_Run_';

for k=0:size(eIdx,2)-1
   fileName = strcat(wDir, filesep,name, string(compose('%04d', k)));
   data = squeeze(double(dicomread(fileName)));
   data = permute(data, [2,3,1]);

   profileData = squeeze(sum(data(integration,integration,:),[1,2]));
   
   [maxV, maxI] = max(profileData);
   [~, r80ind] = min(abs(profileData(maxI:end) - 0.8 * maxV));
   r80ind = r80ind - 1;
   topasR80_air(k+1) = interp1(profileData(maxI + r80ind - 1:maxI + r80ind + 1), ...
             x(maxI + r80ind - 1:maxI + r80ind + 1), 0.8 * maxV);
   
   

%    plot(x, profileData./max(profileData), '.-', 'color','r');
%    hold on;

end
%% Plot topas air vs vacuum
figure;
subplot(2,1,1);
plot([machine.data(eIdx).energy],baseDataRanges(eIdx), 'o-');
hold on;
plot([machine.data(eIdx).energy],topasR80, 'o-');
plot([machine.data(eIdx).energy],topasR80_air, 'o-');

xlabel('Energy [MeV]');
ylabel('Range [mm]');
grid on;
legend('baseData', 'TOPAS Vacuum', 'TOPAS air');
title('E-R relation, generic baseData');
subplot(2,1,2);
plot([machine.data(eIdx).energy],baseDataRanges(eIdx)-topasR80, 'o-');
hold on;
plot([machine.data(eIdx).energy],baseDataRanges(eIdx)-topasR80_air, 'o-');
xlabel('Energy [MeV]');
ylabel('baseData - topas [mm]');
% legend('Vacuum', 'Air');

yyaxis right;
plot([machine.data(eIdx).energy],topasR80-topasR80_air, '--');
ylabel('vacuum -air [mm]');
grid on;
legend('Vacuum', 'Air', 'difference');

%% HIT BaseData
HitBaseData = load(['protons_HITgantry']);

machineHIT = HitBaseData.machine;
baseDataEnergiesHIT = [machineHIT.data(:).energy];

baseDataRangesHITZ = [];
%%%% MCemittance r80 sampling
for k=1:size(baseDataEnergiesHIT,2)
   newDepths = linspace(0,machineHIT.data(k).depths(end),numel(machineHIT.data(k).depths) * 100);
   newDose   = interp1(machineHIT.data(k).depths, machineHIT.data(k).Z, newDepths, 'spline');

   [maxV, maxI] = max(newDose);
   [~, r80ind] = min(abs(newDose(maxI:end) - 0.8 * maxV));
   r80ind = r80ind - 1;
   baseDataRangesHITZ(k) = interp1(newDose(maxI + r80ind - 1:maxI + r80ind + 1), ...
             newDepths(maxI + r80ind - 1:maxI + r80ind + 1), 0.8 * maxV);
end


%% HIT PROTONS: energy-range baseData vs fit
 meanEnergy = @(x) 5.762374661332111e-20 * x^9 - 9.645413625310569e-17 * x^8 + 7.073049219034644e-14 * x^7 ...
                        - 2.992344292008054e-11 * x^6 + 8.104111934547256e-09 * x^5 - 1.477860913846939e-06 * x^4 ...
                        + 1.873625800704108e-04 * x^3 - 1.739424343114980e-02 * x^2 + 1.743224692623838e+00 * x ...
                        + 1.827112816899668e+01;
fittEnergiesHIT = arrayfun(@(r80) meanEnergy(r80),baseDataRangesHITZ);
fittRangeHIT   = interp1(fittEnergiesHIT,baseDataRangesHITZ,baseDataEnergiesHIT);
figure;
subplot(3,1,1);
plot(baseDataEnergiesHIT,baseDataRangesHITZ, 'o');
hold on;
plot(baseDataEnergiesHIT,fittRangeHIT, '.-');
grid on;
xlabel('Energy [MeV]');
ylabel('Range [mm]');
legend('baseData Range', 'fitted Range');
title('proton baseData HIT vs fitted Range');
subplot(3,1,2);
plot(baseDataEnergiesHIT,baseDataRangesHITZ-fittRangeHIT, '.-');
grid on;
xlabel('Energy [MeV]');
ylabel('baseData - fit [mm]');
subplot(3,1,3);
plot(baseDataRangesHITZ,baseDataEnergiesHIT-fittEnergiesHIT, '.-');
grid on;
xlabel('Range [mm]');
ylabel('baseData - fit [MeV]');

%% HIT protons: Z Plots
eIdx = [100, 150, 200];
figure;
for k=1:size(eIdx,2)
   plot(machineHIT.data(eIdx(k)).depths,machineHIT.data(eIdx(k)).Z, '.-');
   xline(baseDataRangesHITZ(eIdx(k)));
   xline(fittRangeHIT(eIdx(k)), 'color', 'r');
   hold on;
end
grid on;
legend('PDD', 'baseData range', 'fit Range');


%% HIT protons: topas simulations

eIdx = [1:20:size(machineHIT.data,2)];
eIdx = eIdx(1:end-2);
energies = baseDataEnergiesHIT(eIdx);
addOffset = true;
wDir = 'D:\matRad_gitHubRemo_MKM_TOPAS\topas\MCrun\BaseData_protons\Results\HIT_MCenergy_TOPAS_vacuum';

name = 'score_BaseData_field1_run1_physicalDose_MCenergy_vacuum_HIT_Run_';

%figure;
x = [1:1500]*0.2 - 0.1;
integration = [40:60];
%Load vaccum profiles
for k=0:size(eIdx,2)-1
   fileName = strcat(wDir, filesep,name, string(compose('%04d', k)));
   data = squeeze(double(dicomread(fileName)));
   data = permute(data, [2,3,1]);

   profileData = squeeze(sum(data(integration,integration,:),[1,2]));
   
   [maxV, maxI] = max(profileData);
   [~, r80ind] = min(abs(profileData(maxI:end) - 0.8 * maxV));
   r80ind = r80ind - 1;
   topasR80_HIT(k+1) = interp1(profileData(maxI + r80ind - 1:maxI + r80ind + 1), ...
             x(maxI + r80ind - 1:maxI + r80ind + 1), 0.8 * maxV);
    if addOffset
      topasR80_HIT(k+1) = topasR80_HIT(k+1) - machineHIT.data(k+1).offset;
   end
   PDD(:,k+1) = profileData;
   %AUC_PDD(k) = sum(profileData).*(x(2) -x(1));
%    plot(x, profileData./max(profileData), '.-', 'color','r');
%    hold on;

end

%Plot Vacuum profiles
figure;
lege = [];
for k=3:size(eIdx,2)-3
   plot(x, PDD(:,k)./max(PDD(:,k)), '.-');
   hold on;
   plot(machineHIT.data(eIdx(k)).depths,machineHIT.data(eIdx(k)).Z./max(machineHIT.data(eIdx(k)).Z), '--', 'color','k');
   xline(baseDataRangesHITZ(eIdx(k)));
   lege = [lege, {['Topas : ', num2str(energies(k))]}, {['baseData : ', num2str(energies(k))]}];
end
xlabel('depth [mm]');
ylabel('dose');
legend(lege);
grid on;
title('Vacuum');

%Load air profiles
wDir = 'D:\matRad_gitHubRemo_MKM_TOPAS\topas\MCrun\BaseData_protons\Results\HIT_MCenergy_TOPAS_air';
name = 'score_BaseData_field1_run1_physicalDose_MCenergy_air_HIT_Run_';

for k=0:size(eIdx,2)-1
   fileName = strcat(wDir, filesep,name, string(compose('%04d', k)));
   data = squeeze(double(dicomread(fileName)));
   data = permute(data, [2,3,1]);

   profileData = squeeze(sum(data(integration,integration,:),[1,2]));
   
   [maxV, maxI] = max(profileData);
   [~, r80ind] = min(abs(profileData(maxI:end) - 0.8 * maxV));
   r80ind = r80ind - 1;
   topasR80_air_HIT(k+1) = interp1(profileData(maxI + r80ind - 1:maxI + r80ind + 1), ...
             x(maxI + r80ind - 1:maxI + r80ind + 1), 0.8 * maxV);
   if addOffset
      topasR80_air_HIT(k+1) = topasR80_air_HIT(k+1) - machineHIT.data(k+1).offset;
   end
end
%% HIT PROTONS: topas plots
figure;
subplot(2,1,1);
plot([machineHIT.data(eIdx).energy],baseDataRangesHITZ(eIdx), 'o-');
hold on;
plot([machineHIT.data(eIdx).energy],topasR80_HIT, 'o-');
plot([machineHIT.data(eIdx).energy],topasR80_air_HIT, 'o-');

xlabel('Energy [MeV]');
ylabel('Range [mm]');
grid on;
legend('baseData', 'TOPAS Vacuum', 'TOPAS air');
title('E-R relation, generic baseData');
subplot(2,1,2);
plot([machineHIT.data(eIdx).energy],baseDataRangesHITZ(eIdx)-topasR80_HIT, 'o-');
hold on;
plot([machineHIT.data(eIdx).energy],baseDataRangesHITZ(eIdx)-topasR80_air_HIT, 'o-');
xlabel('Energy [MeV]');
ylabel('baseData - topas [mm]');
% legend('Vacuum', 'Air');

yyaxis right;
plot([machineHIT.data(eIdx).energy],topasR80_HIT-topasR80_air_HIT, '--');
ylabel('vacuum -air [mm]');
grid on;
legend('Vacuum', 'Air', 'difference');
%% PROTONS MCSquare: load baseData

MCSBaseData = load(['protons_generic_TOPAS_mcSquare.mat']);

MCSmachine = MCSBaseData.machine;
baseDataEnergiesMCS = [MCSmachine.data(:).energy];

%%%% MCemittance r80 sampling
for k=1:size(baseDataEnergiesMCS,2)
   newDepths = linspace(0,MCSmachine.data(k).depths(end),numel(MCSmachine.data(k).depths) * 100);
   newDose   = interp1(MCSmachine.data(k).depths, MCSmachine.data(k).Z, newDepths, 'spline');

   [maxV, maxI] = max(newDose);
   [~, r80ind] = min(abs(newDose(maxI:end) - 0.8 * maxV));
   r80ind = r80ind - 1;
   baseDataRangesMCS(k) = interp1(newDose(maxI + r80ind - 1:maxI + r80ind + 1), ...
             newDepths(maxI + r80ind - 1:maxI + r80ind + 1), 0.8 * maxV);
end
%% PROTONS MCS: fit E-R
 meanEnergy = @(x) 5.762374661332111e-20 * x^9 - 9.645413625310569e-17 * x^8 + 7.073049219034644e-14 * x^7 ...
                        - 2.992344292008054e-11 * x^6 + 8.104111934547256e-09 * x^5 - 1.477860913846939e-06 * x^4 ...
                        + 1.873625800704108e-04 * x^3 - 1.739424343114980e-02 * x^2 + 1.743224692623838e+00 * x ...
                        + 1.827112816899668e+01;
fittEnergiesMCS = arrayfun(@(r80) meanEnergy(r80),baseDataRangesMCS);
fittRangeMCS  = interp1(fittEnergiesMCS,baseDataRangesMCS,baseDataEnergiesMCS);
figure;
subplot(2,1,1);
plot(baseDataEnergiesMCS,baseDataRangesMCS, 'o');
hold on;
plot(baseDataEnergiesMCS,fittRangeMCS, '.-');
grid on;
xlabel('Energy [MeV]');
ylabel('Range [mm]');
legend('baseData Range', 'fitted Range');
title('proton baseData MCS vs fitted Range');
subplot(2,1,2);
plot(baseDataEnergiesMCS,baseDataRangesMCS-fittRangeMCS, '.-');
grid on;
xlabel('Energy [MeV]');
ylabel('baseData - fit [mm]');

%% PROTONS MCS: plot Z profiles
eIdx = [20,30,50,70];
figure;
for k=1:size(eIdx,2)
   plot(MCSmachine.data(eIdx(k)).depths,MCSmachine.data(eIdx(k)).Z, '.-');
   xline(baseDataRangesMCS(eIdx(k)));
   xline(fittRangeMCS(eIdx(k)), 'color', 'r');
   hold on;
end
grid on;
legend('PDD', 'baseData range', 'fit Range');
%% PROTON MCS: TOPAS

eIdx = [1:10:size(MCSmachine.data,2)];
eIdx = eIdx(1:end-2);
energies = baseDataEnergiesMCS(eIdx);
addOffset = false;
wDir = 'D:\matRad_gitHubRemo_MKM_TOPAS\topas\MCrun\BaseData_protons\Results\MCS_MCenergy_TOPAS_vacuum';
name = 'score_BaseData_field1_run1_physicalDose_MCenergy_vacuum_Run_';

%figure;
x = [1:1500]*0.2 - 0.1;
integration = [40:60];
%Load vaccum profiles
PDD = [];
for k=0:size(eIdx,2)-1
   fileName = strcat(wDir, filesep,name, string(compose('%04d', k)));
   data = squeeze(double(dicomread(fileName)));
   data = permute(data, [2,3,1]);

   profileData = squeeze(sum(data(integration,integration,:),[1,2]));
   
   [maxV, maxI] = max(profileData);
   [~, r80ind] = min(abs(profileData(maxI:end) - 0.8 * maxV));
   r80ind = r80ind - 1;
   topasR80_MCS(k+1) = interp1(profileData(maxI + r80ind - 1:maxI + r80ind + 1), ...
             x(maxI + r80ind - 1:maxI + r80ind + 1), 0.8 * maxV);
    if addOffset
      topasR80_MCS(k+1) = topasR80_MCS(k+1) - MCSmachine.data(k+1).offset;
   end
   PDD(:,k+1) = profileData;
   %AUC_PDD(k) = sum(profileData).*(x(2) -x(1));
%    plot(x, profileData./max(profileData), '.-', 'color','r');
%    hold on;

end

%Plot Vacuum profiles
figure;
lege = [];
for k=3:size(eIdx,2)-3
   plot(x, PDD(:,k)./max(PDD(:,k)), '.-');
   hold on;
   plot(MCSmachine.data(eIdx(k)).depths,MCSmachine.data(eIdx(k)).Z./max(MCSmachine.data(eIdx(k)).Z), '--', 'color','k');
   xline(baseDataRangesMCS(eIdx(k)));

   lege = [lege, {['Topas : ', num2str(energies(k))]}, {['baseData : ', num2str(energies(k))]}];
end
xlabel('depth [mm]');
ylabel('dose');
legend(lege);
grid on;

title('Vacuum');

%Load air profiles
wDir = 'D:\matRad_gitHubRemo_MKM_TOPAS\topas\MCrun\BaseData_protons\Results\MCS_MCenergy_TOPAS_air';
name = 'score_BaseData_field1_run1_physicalDose_MCenergy_air_Run_';

for k=0:size(eIdx,2)-1
   fileName = strcat(wDir, filesep,name, string(compose('%04d', k)));
   data = squeeze(double(dicomread(fileName)));
   data = permute(data, [2,3,1]);

   profileData = squeeze(sum(data(integration,integration,:),[1,2]));
   
   [maxV, maxI] = max(profileData);
   [~, r80ind] = min(abs(profileData(maxI:end) - 0.8 * maxV));
   r80ind = r80ind - 1;
   topasR80_air_MCS(k+1) = interp1(profileData(maxI + r80ind - 1:maxI + r80ind + 1), ...
             x(maxI + r80ind - 1:maxI + r80ind + 1), 0.8 * maxV);
   if addOffset
      topasR80_air_MCS(k+1) = topasR80_air_MCS(k+1) - machineMCS.data(k+1).offset;
   end
end

%% PROTONS MCS: Plot topas
figure;
subplot(2,1,1);
plot([MCSmachine.data(eIdx).energy],baseDataRangesMCS(eIdx), 'o-');
hold on;
plot([MCSmachine.data(eIdx).energy],topasR80_MCS, 'o-');
plot([MCSmachine.data(eIdx).energy],topasR80_air_MCS, 'o-');

xlabel('Energy [MeV]');
ylabel('Range [mm]');
grid on;
legend('baseData', 'TOPAS Vacuum', 'TOPAS air');
title('E-R relation, MCS baseData');
subplot(2,1,2);
plot([MCSmachine.data(eIdx).energy],baseDataRangesMCS(eIdx)-topasR80_MCS, 'o-');
hold on;
plot([MCSmachine.data(eIdx).energy],baseDataRangesMCS(eIdx)-topasR80_air_MCS, 'o-');

xlabel('Energy [MeV]');
ylabel('baseData - topas [mm]');
% legend('Vacuum', 'Air');

yyaxis right;
plot([MCSmachine.data(eIdx).energy],topasR80_MCS-topasR80_air_MCS, '--');
ylabel('vacuum -air [mm]');
grid on;
legend('Vacuum', 'Air', 'difference');

%% protons MSC: Topas WATER phantom
eIdx = [1:10:size(MCSmachine.data,2)];
eIdx = eIdx(1:end-2);
energies = baseDataEnergiesMCS(eIdx);
addOffset = false;

wDir = 'D:\matRad_gitHubRemo_MKM_TOPAS\topas\MCrun\BaseData_protons\Results\WATERPhantom\MCS_MCenergy_TOPAS_water';
name = 'score_BaseData_field1_run1_physicalDose_MCenergy_vacuum_WATER_Run_';
%Load vacuum profiles
for k=0:size(eIdx,2)-1
   fileName = strcat(wDir, filesep,name, string(compose('%04d', k)));
   data = squeeze(double(dicomread(fileName)));
   data = permute(data, [2,3,1]);

   profileData = squeeze(sum(data(integration,integration,:),[1,2]));
   
   [maxV, maxI] = max(profileData);
   [~, r80ind] = min(abs(profileData(maxI:end) - 0.8 * maxV));
   r80ind = r80ind - 1;
   topasR80_water_MCS(k+1) = interp1(profileData(maxI + r80ind - 1:maxI + r80ind + 1), ...
             x(maxI + r80ind - 1:maxI + r80ind + 1), 0.8 * maxV);
   if addOffset
      topasR80_water_MCS(k+1) = topasR80_water_MCS(k+1) - machineMCS.data(k+1).offset;
   end
   
   Pdd(:,k+1) = profileData;

end

figure;
lege = [];
for k=3:size(eIdx,2)-3
   plot(x, Pdd(:,k)./max(Pdd(:,k)), '.-');
   hold on;
   plot(MCSmachine.data(eIdx(k)).depths,MCSmachine.data(eIdx(k)).Z./max(MCSmachine.data(eIdx(k)).Z), '--', 'color','k');
   xline(baseDataRangesMCS(eIdx(k)));

   lege = [lege, {['Topas : ', num2str(energies(k))]}, {['baseData : ', num2str(energies(k))]}];
end
xlabel('depth [mm]');
ylabel('dose');

legend(lege);
grid on;
title('Vacuum');

wDir = 'D:\matRad_gitHubRemo_MKM_TOPAS\topas\MCrun\BaseData_protons\Results\WATERPhantom\MCS_MCenergy_TOPAS_air';
name = 'score_BaseData_field1_run1_physicalDose_MCenergy_vacuum_WATER_air_Run_';
% wDir = 'D:\matRad_gitHubRemo_MKM_TOPAS\topas\MCrun\BaseData_protons\Results\WATERPhantom\MCS_MCenergy_TOPAS_airNonCorrected';
% name = 'score_BaseData_field1_run1_physicalDose_MCenergy_airNonCorrected_WATER_Run_';
%Load air profiles
for k=0:size(eIdx,2)-1
   fileName = strcat(wDir, filesep,name, string(compose('%04d', k)));
   data = squeeze(double(dicomread(fileName)));
   data = permute(data, [2,3,1]);

   profileData = squeeze(sum(data(integration,integration,:),[1,2]));
   
   [maxV, maxI] = max(profileData);
   [~, r80ind] = min(abs(profileData(maxI:end) - 0.8 * maxV));
   r80ind = r80ind - 1;
   topasR80_water_MCS_air(k+1) = interp1(profileData(maxI + r80ind - 1:maxI + r80ind + 1), ...
             x(maxI + r80ind - 1:maxI + r80ind + 1), 0.8 * maxV);
   if addOffset
      topasR80_water_MCS_air(k+1) = topasR80_water_MCS_air(k+1) - machineMCS.data(k+1).offset;
   end
end


%% PROTONS MCS: Plot topas WATER PHANTOM
figure;
subplot(2,1,1);
plot([MCSmachine.data(eIdx).energy],baseDataRangesMCS(eIdx), 'o-');
hold on;
plot([MCSmachine.data(eIdx).energy],topasR80_water_MCS, 'o-');
plot([MCSmachine.data(eIdx).energy],topasR80_water_MCS_air, 'o-');

xlabel('Energy [MeV]');
ylabel('Range [mm]');
grid on;
legend('baseData', 'TOPAS Vacuum', 'TOPAS air');
title('E-R relation, MCS baseData Water phantom');
subplot(2,1,2);

plot([MCSmachine.data(eIdx).energy],baseDataRangesMCS(eIdx)-topasR80_water_MCS, 'o-');
hold on;
plot([MCSmachine.data(eIdx).energy],baseDataRangesMCS(eIdx)-topasR80_water_MCS_air, 'o-');

xlabel('Energy [MeV]');
ylabel('baseData - topas [mm]');
% legend('Vacuum', 'Air');

yyaxis right;
plot([MCSmachine.data(eIdx).energy],topasR80_water_MCS-topasR80_water_MCS_air, '--');
ylabel('vacuum -air [mm]');
grid on;
legend('Vacuum', 'Air with dR correction', 'difference');

% wDir = '';
% name = 'score_BaseData_field1_run1_physicalDose_MCenergy_vacuum_WATER_Run_';
% %Load vacuum profiles
% for k=0:size(eIdx,2)-1
%    fileName = strcat(wDir, filesep,name, string(compose('%04d', k)));
%    data = squeeze(double(dicomread(fileName)));
%    data = permute(data, [2,3,1]);
% 
%    profileData = squeeze(sum(data(integration,integration,:),[1,2]));
%    
%    [maxV, maxI] = max(profileData);
%    [~, r80ind] = min(abs(profileData(maxI:end) - 0.8 * maxV));
%    r80ind = r80ind - 1;
%    topasR80_water_MCS(k+1) = interp1(profileData(maxI + r80ind - 1:maxI + r80ind + 1), ...
%              (x(maxI + r80ind - 1:maxI + r80ind + 1)+ 0.1)*2 -0.2, 0.8 * maxV);
%    if addOffset
%       topasR80_water_MCS(k+1) = topasR80_water_MCS(k+1) - machineMCS.data(k+1).offset;
%    end
%    
%    Pdd(:,k+1) = profileData;
% 
% end

%% Carbons
clear all;
meanEnergy = @(x) 11.39 * x^0.628 + 11.24;
ranges = linspace(0,300,1000);
fittedEnergies = arrayfun(@(r80) meanEnergy(r80),ranges);

figure;
plot(fittedEnergies, ranges, '.-');
grid on;
xlabel('E [MeV/u]');
ylabel('Range [mm]');
title('carbon meanEnergy fit');
%% baseData
load('carbon_Generic.mat');
baseDataEnergies = [machine.data(:).energy];
baseDataRanges   = [machine.data(:).range];

%% CARBON: baseData range vs baseData Z r80
baseDataRangesZ = [];
%%%% MCemittance r80 sampling
for k=1:size(baseDataEnergies,2)
   newDepths = linspace(0,machine.data(k).depths(end),numel(machine.data(k).depths) * 100);
   newDose   = interp1(machine.data(k).depths, machine.data(k).Z, newDepths, 'spline');

   [maxV, maxI] = max(newDose);
   [~, r80ind] = min(abs(newDose(maxI:end) - 0.8 * maxV));
   r80ind = r80ind - 1;
   baseDataRangesZ(k) = interp1(newDose(maxI + r80ind - 1:maxI + r80ind + 1), ...
             newDepths(maxI + r80ind - 1:maxI + r80ind + 1), 0.8 * maxV);
end
figure;
subplot(2,1,1);
plot(baseDataEnergies,baseDataRanges, '.-');
hold on;
plot(baseDataEnergies,baseDataRangesZ, '.-');
grid on;
xlabel('Energy [MeV/u]');
ylabel('Range [mm]');
legend('baseData range', 'interpolated r80');
title('carbon baseData vs interp r80');
subplot(2,1,2);
plot(baseDataEnergies,baseDataRanges - baseDataRangesZ, '.-');
grid on;
xlabel('Energy [MeV/u]');
ylabel('baseData - r80 [mm]');

%% CARBON: Z Plots
%%%% Z plot
figure;
plot(machine.data(100).depths, machine.data(100).Z, '.-');
xline(machine.data(100).range);
hold on;
plot(machine.data(10).depths, machine.data(10).Z, '.-');
xline(machine.data(10).range);
grid on;
legend('E = 359.9 MeV/u', 'E = 147.4 [MeV/u]');
title('carbon Z');
xlim([0, 300]);

%% CARBON: baseDataZ range vs fitted range
%%% fitmeanEnergy vs Z baseData
fittedEnergies = arrayfun(@(r80) meanEnergy(r80),baseDataRangesZ);
fittedRanges = interp1(fittedEnergies,baseDataRangesZ,baseDataEnergies, 'spline');

figure;
subplot(2,1,1);
plot(baseDataEnergies,baseDataRangesZ, '.-');
hold on;
plot(baseDataEnergies,fittedRanges, '.-');
grid on;
xlabel('Energy [Mev/u]');
ylabel('Range [mm]');
legend('baseData', 'fitted');
title('carbon baseData range vs fitted range');
subplot(2,1,2);
plot(baseDataEnergies, baseDataRangesZ - fittedRanges, '.-');
grid on;
xlabel('Energy [Mev/u]');
ylabel('baseData - fit [mm]');

%% CARBON: Try fit
F = fittype('a + b*x^c', 'coeff', {'a', 'b', 'c'});
outFit = fit(baseDataRangesZ', baseDataEnergies', F, 'StartPoint', [11.24, 11.39, 0.628]);

%% CARBON: Plot new fit

reFitEnergies = outFit(baseDataRangesZ);
reFitRanges   = interp1(reFitEnergies,baseDataRangesZ,baseDataEnergies);
figure;
subplot(2,1,1);
plot(baseDataEnergies,baseDataRangesZ, '.-');
hold on;
plot(baseDataEnergies, reFitRanges, '.-');
grid on;

xlabel('Energy [MeV/u]');
ylabel('Range [mm]');
subplot(2,1,2);
plot(baseDataEnergies,baseDataRangesZ-reFitRanges, '.-');
grid on;
xlabel('Energy [MeV/u]');
ylabel('Range [mm]');
%% CARBON: baseData range vs fitted range
%This is meaningless, there's no correlation
fittedEnergies_bD = arrayfun(@(r80) meanEnergy(r80),baseDataRanges);
fittedRanges_bD = interp1(fittedEnergies,baseDataRanges,baseDataEnergies, 'spline');

figure;
subplot(2,1,1);
plot(baseDataEnergies,baseDataRanges, '.-');
hold on;
plot(baseDataEnergies,fittedRanges_bD, '.-');
grid on;
xlabel('Energy [Mev/u]');
ylabel('Range [mm]');
legend('baseData', 'fitted');
title('carbon baseData range vs fitted range');
subplot(2,1,2);
plot(baseDataEnergies, baseDataRanges - fittedRanges_bD, '.-');
grid on;
xlabel('Energy [Mev/u]');
ylabel('baseData - fit [mm]');
