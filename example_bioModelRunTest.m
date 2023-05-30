%load generic matRad environment and phantom
matRad_rc;
matRad_cfg = MatRad_Config.instance();
load('BOXPHANTOM.mat');
%% Plan setup

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

%% Instantiate table for output

outcomeStruct = [];
%% Select quantities to be tested
quantitiesOptimization = {'physicalDose', 'RBExD', 'effect'};

%% Select Models to be tested
modelNamesPhotons = {};
modelNamesProtons = {'constRBE', 'MCN', 'WED', 'CAR'};
modelNamesCarbon = {'LEM'};

%% Compute Plan for photons
pln.radiationMode = 'photons';
pln.machine = 'Generic';

stf = matRad_generateStf(ct,cst,pln);
pln.bioParam = [];
for qtOpt = quantitiesOptimization
    %for model = modelNamesPhotons
        currentTable.setupRadiationMode = pln.radiationMode;
        currentTable.setupQuantityOpt   = qtOpt;
        %currentTable.setupModel         = model{1};
        %pln.bioParam = matRad_bioModel(pln.radiationMode,model{1});

        %currentTable.model = pln.bioParam.model;
        
        try
            pln.propOpt.quantityOpt = qtOpt{1};
            dij = matRad_calcPhotonDose(ct,stf,pln,cst,0);
            resultGUI = matRad_fluenceOptimization(dij,cst,pln);
            currentTable.calcBioParametrs = [];
            currentTable.optimized = true;
        catch
            currentTable.optimized = false;
        end
        outcomeStruct = [outcomeStruct, currentTable];
    %end
end



%% Compute Plan for protons
pln.radiationMode = 'protons';
pln.machine = 'Generic';

stf = matRad_generateStf(ct,cst,pln);

for qtOpt = quantitiesOptimization
    for model = modelNamesProtons
        currentTable.setupRadiationMode = pln.radiationMode;
        currentTable.setupQuantityOpt = qtOpt;
        currentTable.setupModel       = model{1};
        

        pln.bioParam = matRad_bioModel(pln.radiationMode,model{1});

        currentTable.model = pln.bioParam.model;
        
        try
            pln.propOpt.quantityOpt = qtOpt{1};
            dij = matRad_calcParticleDose(ct,stf,pln,cst,0);
            resultGUI = matRad_fluenceOptimization(dij,cst,pln);
            currentTable.calcBioParametrs = pln.bioParam.calcBioParameters;
            currentTable.optimized = true;
        catch
            currentTable.optimized = false;
        end
        outcomeStruct = [outcomeStruct, currentTable];
    end
end

% %% Ill combinations
% quantitiesOptimization = {'physicalDose'};
% modelNamesPhotons = {'MCN', 'WED', 'CAR', 'LEM'};
% modelNamesProtons = {'LEM'};
% modelNamesCarbon = {'MCN', 'WED', 'CAR', 'LEM'};
%% Compute Plan for Carbons
pln.radiationMode = 'carbon';
pln.machine = 'Generic';

stf = matRad_generateStf(ct,cst,pln);

for qtOpt = quantitiesOptimization
    for model = modelNamesCarbon
        currentTable.setupRadiationMode = pln.radiationMode;
        currentTable.setupQuantityOpt = qtOpt;
        currentTable.setupModel         = model{1};
        
        pln.bioParam = matRad_bioModel(pln.radiationMode,model{1});

        currentTable.model = pln.bioParam.model;
        
        try
            pln.propOpt.quantityOpt = qtOpt{1};
            dij = matRad_calcParticleDose(ct,stf,pln,cst,0);
            resultGUI = matRad_fluenceOptimization(dij,cst,pln);
            currentTable.calcBioParametrs = pln.bioParam.calcBioParameters;
            currentTable.optimized = true;
        catch
            currentTable.optimized = false;
        end
        outcomeStruct = [outcomeStruct, currentTable];
    end
end

%% Print table for results
outcomeTable = table('Size', [length(outcomeStruct) 3], 'VariableTypes', {'string', 'logical', 'logical'});
outcomeTable.Properties.VariableNames = {'model', 'calcBioParameters', 'optimized'};
%outcomeTable = outcomeStruct;

for k = 1:(length(outcomeStruct))
    rowNames{k} = [outcomeStruct(k).setupRadiationMode, ',', outcomeStruct(k).setupQuantityOpt{1}, ',', outcomeStruct(k).setupModel];
end

outcomeTable.Properties.RowNames = rowNames;
outcomeTable.model               = {outcomeStruct.model}';
outcomeTable.calcBioParameters   = {outcomeStruct.calcBioParametrs}';
outcomeTable.optimized           = {outcomeStruct.optimized}';

%% Rerun model photons
pln.radiationMode = 'photons';
quantityOpt = 'RBExD';
modelName = 'none';
stf = matRad_generateStf(ct,cst,pln);

pln.bioParam = matRad_bioModel(pln.radiationMode, quantityOpt,modelName,pln.machine);
dij = matRad_calcPhotonDose(ct,stf, pln, cst, 0);
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Rerun model protons
pln.radiationMode = 'protons';
quantityOpt = 'RBExD';
modelName = 'constRBE';
stf = matRad_generateStf(ct,cst,pln);

pln.bioParam = matRad_bioModel(pln.radiationMode, quantityOpt,modelName,pln.machine);
dij = matRad_calcParticleDose(ct,stf, pln, cst, 0);
resultGUI = matRad_fluenceOptimization(dij,cst,pln);