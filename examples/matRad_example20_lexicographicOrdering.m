matRad_rc;
matRad_cfg = MatRad_Config.instance();
%% Pln
load('TG119.mat');

% Set plan paramters
pln.radiationMode = 'protons';
pln.machine       = 'Generic';

pln.propDoseCalc.calcLET = 0;
pln.propDoseCalc.engine = 'HongPB';

pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [-45 45];
pln.propStf.couchAngles   = [0, 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

quantityOpt   = 'physicalDose';
modelName     = 'none';


pln.bioParam = matRad_bioModel(pln.radiationMode, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

pln.propDoseCalc.doseGrid.resolution.x = 8; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 8; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 8; % [mm]

%% stf
stf = matRad_generateStf(ct,cst,pln);

%% Dose calc
dij = matRad_calcDoseInfluence(ct,cst,stf,pln);

%% List

% Define the prioritized list of objectives.
% In this case:
%   - (1) Sq. Deviation on the target
%   - (2) EUD on the core
%   - (3) Mean Dose on the body
% Also set a min and max dose constraint for the target structure

% Instantiate the Lexicographic list. This is an object.
PriorityList = matRad_LexicographicList();

% Use the addObjective and addConstraint routines to add the objectives to
% the List. The syntax is the following:
% matRad_LexicographicList.addObjective(priorityIndex, objective, structureIndex, optimizationQuantity, goalValue, robustness);

PriorityList.addObjective(1,DoseObjectives.matRad_SquaredDeviation(100,50),2,'physicalDose',1,'none');
PriorityList.addObjective(2,DoseObjectives.matRad_EUD(100,0),1,'physicalDose',1,'none');
PriorityList.addObjective(3,DoseObjectives.matRad_MeanDose(100,0),3,'physicalDose',1,'none');
PriorityList.addConstraint(DoseConstraints.matRad_MinMaxDose(45,55),2,'physicalDose','none');

% NOTE: The goal value is only used by the 2pec approach. In this case we
% use a nominal lexicographic optimization so the parameter will be
% disregarded. For now, we still need to assign an arbitrary value. For the
% constraint we do not need to set any goal.
%% Optimization
% Define the directory where the intermediate steps for the optimization
% will be saved.
pln.propOpt.saveDir = fullfile(matRad_cfg.primaryUserFolder, 'TG119_testLexicographicOptimization');

% Define slack quantity. This relaxes the constraint posed on the previous
% objectives.
pln.propOpt.slack = 1.10;

% Increase max iterations so no plan stops prematurely
matRad_cfg.defaults.propOpt.maxIter = 10000;


% Run!
[optimizedList, resultGUIs] = matRad_lexicographicOptimization(pln,dij,cst,PriorityList);

%% Visualize

steps = fieldnames(resultGUIs);
for i=1:numel(steps)
    curCube = resultGUIs.(steps{i}).physicalDose;
    currDVH = matRad_calcDVH(cst,curCube);
    f = figure('WindowState','maximized');
    %f.Position(1:2) = f.Position(1:2)-300;

    nexttile();
    matRad_plotSlice(ct,'dose', curCube, 'axesHandle', gca(), 'cst', cst,'cubeIdx',1, 'plane',3, 'slice', 70, 'doseWindow',[0,2]);
    title(sprintf('Optimization step: %d', i));
    

    nexttile();
    for j=1:size(cst,1)
        plot(currDVH(j).doseGrid, currDVH(j).volumePoints, '-', 'LineWidth',2.5,'DisplayName',cst{j,3});
        hold on;
        xlim([0,2]);
    end
    grid on;
    xlabel('Dose [Gy]');
    ylabel('Volume [%]');
    legend();
end
