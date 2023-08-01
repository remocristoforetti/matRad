matRad_rc;
matRad_cfg = MatRad_Config.instance();
matRad_cfg.propOpt.defaultMaxIter = 50000;
load 'TG119.mat'
cst{1,6}{1}.parameters = {20};

cst{1,6}{1}.penalty = 500;
cst{3,6}{1}.parameters = {20};
cst{1,6}{1}.penalty = 300;
% 
% meta information for treatment plan (1) 
pln(1).numOfFractions  = 15;
pln(1).radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln(1).machine         = 'Generic';

% beam geometry settings
pln(1).propStf.bixelWidth      = 3; % [mm] / also corresponds to lateral spot spacing for particles
pln(1).propStf.gantryAngles    = [45 315]; % [?] ;
pln(1).propStf.couchAngles     = zeros(numel(pln(1).propStf.gantryAngles),1); % [?] ; 
pln(1).propStf.numOfBeams      = numel(pln(1).propStf.gantryAngles);
pln(1).propStf.isoCenter       = ones(pln(1).propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln(1).propDoseCalc.calcLET = 0;

pln(1).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln(1).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln(1).propOpt.spatioTemp      = 0;
pln(1).propOpt.STscenarios     = 2;
%pln(1).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln(1).propDoseCalc.doseGrid.resolution.x = 8; % [mm]
pln(1).propDoseCalc.doseGrid.resolution.y = 8; % [mm]
pln(1).propDoseCalc.doseGrid.resolution.z = 8; % [mm]
% pln(1).propDoseCalc.doseGrid.resolution = ct.resolution;
quantityOpt  = 'physicalDose';     % options: physicalDose, effect, RBExD
%=======================================> Model check error in bioModel
modelName    = 'none';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln(1).bioParam = matRad_bioModel(pln(1).radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln(1).multScen = matRad_multScen(ct,scenGenType);
%pln(1).multScen.includeNominalScenario = 1;
%pln(1).multScen.nSamples = 3;
% 
% meta information for treatment plan (2) 
pln(2).numOfFractions  = 15;
pln(2).radiationMode   = 'photons';           % either photons / protons / helium / carbon
pln(2).machine         = 'Generic';

% beam geometry settings
pln(2).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln(2).propStf.gantryAngles    = [0:50:359]; % [?] ;
pln(2).propStf.couchAngles     = zeros(numel(pln(2).propStf.gantryAngles),1);  % [?] ; 
pln(2).propStf.numOfBeams      = numel(pln(2).propStf.gantryAngles);
pln(2).propStf.isoCenter       = ones(pln(2).propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln(2).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln(2).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln(2).propOpt.spatioTemp      = 0;
pln(2).propOpt.STscenarios     = 5;
%pln(2).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln(2).propDoseCalc.doseGrid.resolution.x = 8; % [mm]
pln(2).propDoseCalc.doseGrid.resolution.y = 8; % [mm]
pln(2).propDoseCalc.doseGrid.resolution.z = 8; % [mm]
% pln(2).propDoseCalc.doseGrid.resolution = ct.resolution;

quantityOpt  = 'physicalDose';     % options: physicalDose, effect, RBExD
modelName    = 'none';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln(2).bioParam = matRad_bioModel(pln(2).radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln(2).multScen = matRad_multScen(ct,scenGenType);
%pln(2).multScen.includeNominalScenario = 1;
%pln(2).multScen.nSamples = 2;

% prepping cst 
% placing alpha/beta ratios in cst{:,6},
% different alpha beta ration for each obj of a structure  
sparecst = 0;

cst = matRad_prepCst(cst, sparecst);


%% Plan Wrapper
plnJO = matRad_plnWrapper(pln,'nominal');
%TODO: properly handle FLAGS for robustness in fluenceOpt
%% Stf Wrapper
stf = matRad_stfWrapper(ct,cst,plnJO);
%% Dij Calculation
dij = matRad_calcCombiDose(ct,stf,plnJO,cst,false);
%% Fluence optimization
% for voiIdx=1:size(cst,1)
% 
%     for objIdx=1:size(cst{voiIdx,6},1)
%         cst{voiIdx,6}{objIdx}.robustness = 'STOCH'; 
%     end
% end
resultGUI_nominal = matRad_fluenceOptimizationJO(dij,cst,plnJO);

%% Robust
scenGenType  = 'rndScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

pln(1).multScen = matRad_multScen(ct,scenGenType);
pln(1).multScen.includeNominalScenario = 1;

pln(1).multScen.nSamples = 7;

scenGenType  = 'rndScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          
pln(2).multScen = matRad_multScen(ct,scenGenType);
pln(2).multScen.includeNominalScenario = 1;
pln(2).multScen.nSamples = 5;


plnJO = matRad_plnWrapper(pln,'allVsall');

% Stf Wrapper
stf = matRad_stfWrapper(ct,cst,plnJO);
% Dij Calculation

dij = matRad_calcCombiDose(ct,stf,plnJO,cst,false);
%% Fluence optimization
for voiIdx=1:size(cst,1)

    for objIdx=1:size(cst{voiIdx,6},1)
        cst{voiIdx,6}{objIdx}.robustness = 'STOCH'; 
    end
end
resultGUI = matRad_fluenceOptimizationJO(dij,cst,plnJO);

%% Visualization
slice = 59;
qtOpt = plnJO.bioParam.quantityOpt;
modalities = {'proton', 'photon'};
plan_JO = pln(1).numOfFractions * resultGUI{1}.(qtOpt) + pln(2).numOfFractions * resultGUI{2}.(qtOpt);

if any(plnJO.propOpt.spatioTemp)
    f = figure;
    protonScenarios = plnJO.propOpt.STscenarios(1);
    if protonScenarios>1
        for scenarioIdx = 1:protonScenarios
            subplot(3,max([plnJO.propOpt.STscenarios])+1,scenarioIdx);
            imagesc(resultGUI{1}.([qtOpt, '_STscenario_', num2str(scenarioIdx)])(:,:,slice));

            matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);
            title(['proton plan scenario ', num2str(scenarioIdx)]);
            colorbar();
        end
    end
    subplot(3,max([plnJO.propOpt.STscenarios])+1,max([plnJO.propOpt.STscenarios])+1);
    imagesc(resultGUI{1}.(qtOpt)(:,:,slice));
    matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);
    title('total proton plan');
    colorbar();
   
    photonScenarios = plnJO.propOpt.STscenarios(2);
    if photonScenarios>1
        for scenarioIdx = 1:photonScenarios
            subplot(3,max([plnJO.propOpt.STscenarios])+1,max([plnJO.propOpt.STscenarios])+1 + scenarioIdx);
            imagesc(resultGUI{2}.([qtOpt, '_STscenario_', num2str(scenarioIdx)])(:,:,slice));
            matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);
            title(['photon plan scenario ', num2str(scenarioIdx)]);
            colorbar();
        end
    end
    subplot(3,max([plnJO.propOpt.STscenarios])+1,2*(max([plnJO.propOpt.STscenarios])+1));
    imagesc(resultGUI{2}.(qtOpt)(:,:,slice));
    matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);
    title('total photon plan');
    colorbar();
    subplot(3,max([plnJO.propOpt.STscenarios])+1,3*(max([plnJO.propOpt.STscenarios])+1));
    imagesc(plan_JO(:,:,slice));
    matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);
    title('total plan');
    colorbar();
else
%%% Plot photon-proton single scenarios total plan
    f = figure;

    for modalityIdx =1:plnJO.numOfModalities
        subplot(1,3,modalityIdx);
        imagesc(resultGUI{modalityIdx}.(qtOpt)(:,:,slice));
        matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);
        title([ modalities{modalityIdx}, ' plan']);
    end

    subplot(1,3,3);
    imagesc(plan_JO(:,:,slice));
    matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);

    title('Total plan');
end
%% Visualization robust
slice = 59;
quantityOpt = 'physicalDose';

if exist('resultGUI_nominal', 'var')

    nominalPlan_photon = resultGUI_nominal{2}.(quantityOpt).*pln(2).numOfFractions;
    nominalPlan_proton = resultGUI_nominal{1}.(quantityOpt).*pln(1).numOfFractions;
    nominalPlan        = pln(1).numOfFractions * resultGUI_nominal{1}.(quantityOpt) + pln(2).numOfFractions * resultGUI_nominal{2}.(quantityOpt);

    Plan_photon = resultGUI{2}.(quantityOpt).*pln(2).numOfFractions;
    Plan_proton = resultGUI{1}.(quantityOpt).*pln(1).numOfFractions;
    Plan        = pln(1).numOfFractions * resultGUI{1}.(quantityOpt) + pln(2).numOfFractions * resultGUI{2}.(quantityOpt);

    f = figure;
    subplot(2,3,1);
    imagesc(nominalPlan_proton(:,:,slice));
    matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);
    title('proton plan nominal');
    colorbar();

    subplot(2,3,2);
    imagesc(nominalPlan_photon(:,:,slice));
    matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);
    title(['photon plan nominal']);
    colorbar();

    subplot(2,3,3);
    imagesc(nominalPlan(:,:,slice));
    matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);
    title(['Nominal plan']);

    colorbar();

    subplot(2,3,4);
    imagesc(Plan_proton(:,:,slice));
    matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);
    title(['proton robust plan']);
    colorbar();

    subplot(2,3,5);
    imagesc(Plan_photon(:,:,slice));
    matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);

    title(['photon robust plan']);
    colorbar();

    subplot(2,3,6);
    imagesc(Plan(:,:,slice));
    matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);
    title(['Robust plan']);

    colorbar();

else

end