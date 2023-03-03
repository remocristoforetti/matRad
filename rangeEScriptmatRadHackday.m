%%% Test E-R fit




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
