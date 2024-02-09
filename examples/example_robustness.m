matRad_rc;

load('TG119.mat'); 

%% Pln

pln.radiationMode = 'protons';
pln.machine       = 'generic';
pln.propDoseCalc.calcLET = 0;

pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [0];
pln.propStf.couchAngles   = zeros(size(pln.propStf.gantryAngles));
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;

pln.propSeq.runSequencing = 0;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 8; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 8; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 8; % [mm]

pln.multScen = matRad_multScen(ct, 'rndScen');
pln.bioParam = matRad_bioModel(pln.radiationMode, 'physicalDose', 'none');
%% Stf

stf = matRad_generateStf(ct,cst,pln);

%% analytical dij

dij = matRad_calcParticleDose(ct,stf,pln,cst);


%% cst
for i=1:size(cst,1)
    for j=1:numel(cst{i,6})
        cst{i,6}{j}.robustness = 'OWC';
    end
end

%% Fluence opt
resultGUI = matRad_fluenceOptimization(dij,cst,pln);