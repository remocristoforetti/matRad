%% Script 1: Build an elongated box phantom for single-PB / small-SOBP demos

%% set matRad runtime configuration
clear;
matRad_rc; %If this throws an error, run it from the matRad root directory first to set the paths

matRad_cfg = MatRad_Config.instance();

%% Phantom / target geometry
ctDim        = [200 70 70];    % [x y z] number of voxels
ctResolution = [1 1 1];        % [x y z] mm per voxel -> phantom spans 200 x 50 x 50 mm^3
numOfCtScen  = 1;              % single (static) ct scenario, no motion

phantomHU    = 0;

targetDim    = [50 50 50];     % [x y z] mm, small cubic target
targetOffset = [50 0 0];       % [x y z] mm, shifted from the phantom center along x
                               % (beam will later enter from the -x face and travel in +

targetHU     = phantomHU;
builder = matRad_PhantomBuilder(ctDim,ctResolution,numOfCtScen);

%% Define the dose objective for the target

targetObjective = struct(DoseObjectives.matRad_SquaredDeviation(1,60));

%% Add the volumes of interest

builder.addBoxTarget('TARGET',targetDim,'HU',targetHU,'offset',targetOffset,'objectives',targetObjective);
builder.addBoxOAR('PHANTOM',ctDim,'HU',phantomHU);

[ct,cst] = builder.getctcst();
ct = matRad_getWorldAxes(ct);

%% Sanity-check figure: visualize the phantom & target before any dose calc
isoCenter = matRad_getIsoCenter(cst,ct,0);
sliceIdx  = matRad_world2cubeIndex(isoCenter,ct);

figure('Name','Box phantom','Color','w');
matRad_plotSlice(ct,'axesHandle',gca,'cst',cst,'plane',3,'slice',sliceIdx(3), ...
    'boolPlotLegend',true, 'title', ' ');

%% Save for reuse in the following scripts (single PB / SOBP demos)
outFile = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','phantomBuilder','phantom_SOBP_2.mat');
save(outFile,'ct','cst');
