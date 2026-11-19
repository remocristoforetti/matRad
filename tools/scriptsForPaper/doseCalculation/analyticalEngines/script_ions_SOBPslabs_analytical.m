%% Script: Hong dose calculation and SOBP optimization on the slab phantom for helium, carbon and oxygen
clear;
matRad_rc;
matRad_cfg = MatRad_Config.instance();

ions = {'helium','carbon','oxygen'};

outDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','doseCalculation','analyticalEngines');
if ~isfolder(outDir)
    mkdir(outDir);
end

for i = 1:numel(ions)
    load(fullfile(matRad_cfg.primaryUserFolder,'dataPaper','phantomBuilder','phantom_SOBP_slabs_2.mat'),'ct','cst');

    cst{strcmp(cst(:,2),'TARGET'),6}{1}  = DoseObjectives.matRad_SquaredDeviation(100,60);
    cst{strcmp(cst(:,2),'PHANTOM'),6}{1} = DoseObjectives.matRad_SquaredOverdosing(100,50);

    clear pln;
    pln.radiationMode   = ions{i};
    pln.machine         = 'Generic';
    pln.bioModel        = 'none';
    pln.multScen        = 'nomScen';
    pln.numOfFractions  = 30;

    pln.propStf.gantryAngles = 270;
    pln.propStf.couchAngles  = 0;
    pln.propStf.bixelWidth   = 3;
    pln.propStf.isoCenter    = matRad_getIsoCenter(cst,ct,0);

    pln.propDoseCalc.engine = 'HongPB';
    pln.propDoseCalc.doseGrid.resolution.x = 1.5;
    pln.propDoseCalc.doseGrid.resolution.y = 1.5;
    pln.propDoseCalc.doseGrid.resolution.z = 1.5;

    pln.propOpt.runDAO = false;

    stf = matRad_generateStf(ct,cst,pln);
    dij = matRad_calcDoseInfluence(ct,cst,stf,pln);
    resultGUI = matRad_fluenceOptimization(dij,cst,pln);

    doseSOBP = resultGUI.physicalDose;
    w        = resultGUI.w;

    save(fullfile(outDir,[ions{i} '_SOBPslabs_analytical_2.mat']),'ct','cst','pln','stf','w','doseSOBP');
end
