%% generate fake matRad information
matRad_rc;

load('BOXPHANTOM.mat');

pln.radiationMode = 'protons';

pln.machine       = 'generic_MCsquare';

pln.propDoseCalc.calcLET = 0;

pln.numOfFractions        = 10;
pln.propStf.gantryAngles  = [0];
pln.propStf.couchAngles   = zeros(size(pln.propStf.gantryAngles));
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 1; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 1; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 1; % [mm]


pln.bioParam = matRad_bioModel(pln.radiationMode, 'physicalDose','none');

pln.multScen = matRad_multScen(ct, 'nomScen');
%% Stf
%stf = matRad_generateStf(ct,cst,pln);

load([pln.radiationMode,'_', pln.machine]);
x = 0;%[-100:30:100];%[-100:30:100];
y = 0;%[-100:30:100];%[-100:30:100];
stf = matRad_generateStfSpotGridForTesting(ct,cst,pln,machine.data(57).energy, x,y);

weights = linspace(1,1,stf.totalNumOfBixels)';
resultGUI.w = weights;
%% FRED Calculation
matRad_cfg = MatRad_Config.instance();
matRad_cfg.propDoseCalc.defaultNumHistoriesPerBeamlet = 100000*sum([stf.totalNumOfBixels]);%1000*dij_PB.totalNumOfBixels;%floor((1000*dij_PB.totalNumOfBixels)/sum(resultGUI.w));

pln.propDoseCalc.engine = 'FRED';

pln.propDoseCalc.exportCalculation = true;
pln.propDoseCalc.sourceModel = 'emittance';

resultGUI_FRED = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);

%% TOPAS
% pln.propDoseCalc.engine = 'TOPAS';
% resultGUI_TOPAS = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);


%% Read manual simulation
sDir = 'C:\r408i_data\r408i_data\matRad_dev\FRED\MCrun_singleSpot_MCsquareGeneric\out\reg\';
doses = [];

dist = [0.2,25,50,75,100,125,150,175,200];
for distIdx=1:numel(dist)

    pName = [sDir, 'Phantom_', num2str(dist(distIdx)), '\'];
    doses(:,:,distIdx) = matRad_readMhd(pName, 'Dose.mhd');

end


figure;

%tiledlayout(1,numel());
for i=1:numel(dist)
    nexttile;


    imagesc(doses(:,:,i));

end



%% Profiles
figure;
x = linspace(-250,250,1000);
leg = [];

for i=1:numel(dist)
    currProf = squeeze(sum(doses(500,:,i), 1));
    plot(x,currProf, '.-');
    hold on;
    leg = [leg, {num2str(dist(i))}];
end

grid on;
grid minor;
xlim([-100,100]);
legend(leg);
xlabel('mm', 'FontSize',17);
ylabel('Dose', 'FontSize',17);

%% Fits
figure;
leg = [];
f=[];
for i=1:numel(dist)
   currProf = squeeze(sum(doses(:,:,i), 1));

   f = fit(x',currProf','gauss1');
   sigm(i) = f.c1/sqrt(2);
   plot(x,currProf, '.');
   hold on;

   p = plot(x,f(x),'-');
   set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
   leg = [leg, {num2str(dist(i))}];
end

grid on;
grid minor;
xlim([-100,100]);

legend(leg);
xlabel('mm', 'FontSize',17);
ylabel('Dose', 'FontSize',17);

%% plot

eIdx = 57;
figure;
tiledlayout(2,1);
nexttile;
plot(machine.data(eIdx).initFocus.dist -(machine.meta.SAD - machine.meta.BAMStoIsoDist), machine.data(eIdx).initFocus.sigma, '.-');
%plot(machine.data(37).initFocus.dist -(machine.meta.SAD), machine.data(37).initFocus.sigma, '.-');
hold on;
plot(dist*10, sigm, '.-');
grid on;
legend('machine', 'FRED');
xlim([0, 2000]);
ylabel('sigma [mm]', 'FontSize',17);
x_diff = linspace(min(machine.data(eIdx).initFocus.dist -(machine.meta.SAD - machine.meta.BAMStoIsoDist)), max(machine.data(eIdx).initFocus.dist -(machine.meta.SAD - machine.meta.BAMStoIsoDist)), 100);
machine_diff = interp1(machine.data(eIdx).initFocus.dist -(machine.meta.SAD - machine.meta.BAMStoIsoDist), machine.data(eIdx).initFocus.sigma,x_diff, 'spline');
fred_diff = interp1(dist*10, sigm, x_diff, 'spline');

nexttile;
plot(x_diff, 100*(machine_diff - fred_diff)./machine_diff, '.-');

xlim([0, 2000]);
ylim([-30, 2.5]);
grid on;
grid  minor;

xlabel('Distance from BAMs [mm]', 'FontSize',17);
ylabel('Residues [%]');