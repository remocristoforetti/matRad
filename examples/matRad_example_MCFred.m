matRad_rc;
load('TG119.mat');

pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

pln.propOpt.bioOptimization = 'physicalDose';

pln.propDoseCalc.engine = 'FRED';

pln.propDoseCalc.calcLET = 0;


pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [0];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 3;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

testStfSingleField;
%stf = matRad_generateStf(ct,cst,pln);

% %% add random weights
% for i=1:numel(stf)
%     for j=1:stf(i).numOfRays
%         stf(i).ray(j).weight = 3*rand(1,stf(i).numOfBixelsPerRay(j),'double');
%     end
% end

%% Dose Calculation

dij = matRad_calcDoseInfluence(ct,cst,stf,pln);

%% read

cube = matRad_readMhd('FRED/MCrun/out', 'Dose.mhd');


%% Visualize
figure;
subplot(1,3,1);
matRad_plotSliceWrapper(gca(),ct,cst,1,cube,1,1);

subplot(1,3,2);
matRad_plotSliceWrapper(gca(),ct,cst,1,cube,2,83);

subplot(1,3,3);
matRad_plotSliceWrapper(gca(),ct,cst,1,cube,3,83);