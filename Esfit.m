BP = dicomread('C:\r408i_data\r408i_data\BaseData_carbon_spectra\SingleEnergyNoEs\Output\PDD\Phantom\physicalDose_Run_0000.dcm');
BP = double(squeeze(BP));
BP = smooth(squeeze(sum(BP, [2,3])));

depthResolution = 0.1;

x = [0:size(BP,1)-1]*depthResolution + depthResolution/2;
% figure;
% plot(x, BP, '.-');
%% Cut central part

window =20; % mm

[~,peakPosIdx] = max(BP);
peakPos = x(peakPosIdx);

newDepthsIdx(1) = find(x<peakPos - window/2, 1, 'last');

newDepthsIdx(2) = find(x<peakPos + window/2, 1, 'last');
cutDepths = x([newDepthsIdx(1):newDepthsIdx(2)]);
cutPDD    = BP([newDepthsIdx(1):newDepthsIdx(2)]);

figure;
plot(cutDepths, cutPDD, 'o-');
%% Oversample

newDepths = linspace(cutDepths(1), cutDepths(end), 10*size(cutDepths,2));
newPDD    = interp1(cutDepths, cutPDD, newDepths, 'spline');

hold on;
plot(newDepths, newPDD, '.-');

%% Get range

[maxV, maxI] = max(newPDD);

[~, r80ind] = min(abs(newPDD(maxI:end) - 0.8 * maxV));
r80ind = r80ind - 1;
r80 = interp1(newPDD(maxI + r80ind - 1:maxI + r80ind + 1),newDepths(maxI + r80ind - 1:maxI + r80ind + 1), 0.8 * maxV);
meanEnergy = @(x) 11.39 * x^0.628 + 11.24;

%% Generate shifts
shifts = 0.5*[-10:7]; % in mm
shiftedR80 = r80 + shifts; %mm


Ebins = arrayfun(@(r) meanEnergy(r), shiftedR80); % these are topas energies
Erange = (Ebins(end)-Ebins(1))*12;

%% Get machine data
load('carbon_Generic.mat');

eIdx = 31;
figure;
for k =1:size(shiftedR80,2)

    plot(newDepths + shifts(k), newPDD./max(newPDD), '.-', 'color', 'k');

    hold on;
end
plot([machine.data(eIdx).depths], machine.data(eIdx).Z./max(machine.data(eIdx).Z), '.-', 'color', 'r');

%% Gen data for minimization
optDepths = newDepths;


PDDs = zeros(size(shifts,2),size(optDepths,2));
for k=1:size(shifts,2)

    PDDs(k,:) = interp1(newDepths+shifts(k), newPDD./max(newPDD), optDepths, 'spline');
end

calcLinComb = @(w,bPs) w*bPs;

w0 = ones(1,size(shifts,2));
objective = interp1(machine.data(eIdx).depths, machine.data(eIdx).Z./max(machine.data(eIdx).Z),optDepths,'spline');

resW = lsqcurvefit(@(w,D) calcLinComb(w,D), w0, PDDs,objective, zeros(size(w0)));

resultPDD = calcLinComb(resW,PDDs);

figure;
plot(optDepths,resultPDD,'.-');
hold on;
plot(optDepths, objective, '--');
for k = 1:size(shifts,2)
    plot(optDepths, resW(k).*PDDs(k,:), '--');
end

%% plot weights

figure;
plot(Ebins, resW./sum(resW), '.-');
%% Gauss fit on weights
resF = fit(Ebins', resW', 'gauss1');

sigma = sqrt(resF.c1/2);

energySpread = (sigma*100)/(resF.b1);