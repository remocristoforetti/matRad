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
%% PROTON MCS: TOPAS Schneider phantom

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

%% PROTONS MCS: Plot topas schneider phantom
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
integration = [40:60];
x = [1:1500]*0.2 - 0.1;

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

% wDir = 'D:\matRad_gitHubRemo_MKM_TOPAS\topas\MCrun\BaseData_protons\Results\WATERPhantom\MCS_MCenergy_TOPAS_air';
% name = 'score_BaseData_field1_run1_physicalDose_MCenergy_vacuum_WATER_air_Run_';
wDir = 'D:\matRad_gitHubRemo_MKM_TOPAS\topas\MCrun\BaseData_protons\Results\WATERPhantom\MCS_MCenergy_TOPAS_airNonCorrected';
name = 'score_BaseData_field1_run1_physicalDose_MCenergy_airNonCorrected_WATER_Run_';
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
legend('Vacuum', 'Air', 'difference');

%% PROTONS MCS: WATER Phantom dR correction
wDir = 'D:\matRad_gitHubRemo_MKM_TOPAS\topas\MCrun\BaseData_protons\Results\WATERPhantom\MCS_MCenergy_TOPAS_air';
name = 'score_BaseData_field1_run1_physicalDose_MCenergy_vacuum_WATER_air_Run_';
%Load vacuum profiles
for k=0:size(eIdx,2)-1
   fileName = strcat(wDir, filesep,name, string(compose('%04d', k)));
   data = squeeze(double(dicomread(fileName)));
   data = permute(data, [2,3,1]);

   profileData = squeeze(sum(data(integration,integration,:),[1,2]));
   
   [maxV, maxI] = max(profileData);
   [~, r80ind] = min(abs(profileData(maxI:end) - 0.8 * maxV));
   r80ind = r80ind - 1;
   topasR80_water_MCS_airCorrected(k+1) = interp1(profileData(maxI + r80ind - 1:maxI + r80ind + 1), ...
             x(maxI + r80ind - 1:maxI + r80ind + 1), 0.8 * maxV);
   if addOffset
      topasR80_water_MCS_airCorrected(k+1) = topasR80_water_MCS_airCorrected(k+1) - machineMCS.data(k+1).offset;
   end
end

%% PROTONS MCS: WATER PHANTOM dR Correction, plot
figure;

subplot(2,1,1);
plot([MCSmachine.data(eIdx).energy],baseDataRangesMCS(eIdx), 'o-');
hold on;
plot([MCSmachine.data(eIdx).energy],topasR80_water_MCS, 'o-');
plot([MCSmachine.data(eIdx).energy],topasR80_water_MCS_air, 'o-');
plot([MCSmachine.data(eIdx).energy],topasR80_water_MCS_airCorrected, 'o-');
xlabel('Energy [MeV]');
ylabel('Range [mm]');
grid on;
legend('baseData', 'TOPAS Vacuum', 'TOPAS air', 'TOPAS air dR correction');
title('E-R relation, MCS baseData Water phantom');
subplot(2,1,2);

plot([MCSmachine.data(eIdx).energy],baseDataRangesMCS(eIdx)-topasR80_water_MCS, 'o-');
hold on;
plot([MCSmachine.data(eIdx).energy],baseDataRangesMCS(eIdx)-topasR80_water_MCS_air, 'o-');
plot([MCSmachine.data(eIdx).energy],baseDataRangesMCS(eIdx)-topasR80_water_MCS_airCorrected, 'o-');

xlabel('Energy [MeV]');
ylabel('baseData - topas [mm]');
% legend('Vacuum', 'Air');

yyaxis right;
plot([MCSmachine.data(eIdx).energy],topasR80_water_MCS-topasR80_water_MCS_airCorrected, '--');
ylabel('vacuum - air (dR correction) [mm]');
grid on;
legend('Vacuum', 'Air','dR Correction' ,'difference');
