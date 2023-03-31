%% Get simulation files
[filenames, dirPath] = uigetfile('.dcm', 'MultiSelect', 'on');

nEnergies = size(filenames,2);
load('carbon_HITgantry.mat');

eIdx = [1:nEnergies];
initialGuess = 1.5;
%% do
for k= 1:nEnergies
   [ES(k),E(k), costFunctionValue(k), exitFlag(k)] = getEnergySpread([dirPath,filenames{k}],machine,eIdx(k),30,initialGuess, 0);
   initialGuess = (initialGuess + ES(k))/2;
   k

end



%%  plot
figure;
plot([machine.data(eIdx(1:k-1)).energy], ES(1:k-1), '.-');

figure;
plot([machine.data(eIdx(1:k-1)).energy], costFunctionValue(1:k-1), '.-');


figure;
plot([machine.data(eIdx(1:k-1)).energy], exitFlag(1:k-1), '.-');

%% ES LUT

EStable.meta.description = 'Energy spread LUT computed for HITgantry base data. 31-3-2023';
EStable.data.E = E;
EStable.data.energySpread = ES;


save('ES.mat', 'EStable');

%% Read out Topas data
[filenames, dirPath] = uigetfile('.dcm', 'MultiSelect', 'on');
%% Load
for k=1:size(filenames,2)
    BP = dicomread([dirPath,filenames{k}]);
    BP = squeeze(double(BP));
    BPs{k} = squeeze(sum(BP,[2,3]));
end
x = [0:size(BPs{1},1)-1]*0.1 + 0.05;
load('carbon_HITgantry.mat');

%% Plot
simEnergies = [machine.data([1:20:end]).energy];
[~,eIdx] = intersect([machine.data.energy], simEnergies);
k=4;



figure;
plot(machine.data(eIdx(k)).depths,machine.data(eIdx(k)).Z./max(machine.data(eIdx(k)).Z), '.-');
hold on;
plot(x, BPs{k}./max(BPs{k}), '.-');
grid on;
grid minor;