%% Script 1: Build the SOBP box phantom with heterogeneous slabs in the entrance channel

%% set matRad runtime configuration
clear;
matRad_rc; %If this throws an error, run it from the matRad root directory first to set the paths
matRad_cfg = MatRad_Config.instance();

%% Phantom / target geometry (identical to matRad_script01_buildPhantom_SOBP.m)
ctDim        = [200 70 70];    % [x y z] number of voxels
ctResolution = [1 1 1];        % [x y z] mm per voxel -> phantom spans 200 x 50 x 50 mm^3
numOfCtScen  = 1;              % single (static) ct scenario, no motion

phantomHU    = 0;

targetDim    = [50 50 50];     % [x y z] mm
targetOffset = [50 0 0];       % [x y z] mm, target spans x = 125..175 mm (beam enters at x = 0)
targetHU     = phantomHU;

%% Slabs in the entrance channel
% Each slab spans the full y/z cross-section of the phantom and sits
% between the entrance face and the target.
slabThickness = 15;                         % [mm] along the beam direction (x)
slabDim       = [slabThickness ctDim(2) ctDim(3)];

slabLowHU     = -500;                       % lower density than water
slabLowOffset = [-60 0 0];                  % spans x = 32.5..47.5 mm

slabHighHU     = 500;                       % higher density than water
slabHighOffset = [-20 0 0];                 % spans x = 72.5..87.5 mm

builder = matRad_PhantomBuilder(ctDim,ctResolution,numOfCtScen);

%% Define the dose objective for the target
targetObjective = struct(DoseObjectives.matRad_SquaredDeviation(1,60));

%% Add the volumes of interest
% getctcst() assigns HU in reverse definition order, so the first-added
% VOI wins overlaps: the slabs must be added before the full-cube PHANTOM,
% otherwise its HU = 0 would overwrite them.
builder.addBoxTarget('TARGET',targetDim,'HU',targetHU,'offset',targetOffset,'objectives',targetObjective);
builder.addBoxOAR('SlabLow',slabDim,'HU',slabLowHU,'offset',slabLowOffset);
builder.addBoxOAR('SlabHigh',slabDim,'HU',slabHighHU,'offset',slabHighOffset);
builder.addBoxOAR('PHANTOM',ctDim,'HU',phantomHU);

[ct,cst] = builder.getctcst();
ct = matRad_getWorldAxes(ct);

%% Sanity-check figure: visualize the phantom & target before any dose calc
isoCenter = matRad_getIsoCenter(cst,ct,0);
sliceIdx  = matRad_world2cubeIndex(isoCenter,ct);

figure('Name','Box phantom with slabs','Color','w');
matRad_plotSlice(ct,'axesHandle',gca,'cst',cst,'plane',3,'slice',sliceIdx(3), ...
    'boolPlotLegend',true, 'title', ' ');

figure('Name','HU along the central beam axis','Color','w');
plot(ct.x,squeeze(ct.cubeHU{1}(sliceIdx(1),:,sliceIdx(3))),'LineWidth',2);
grid on; box on;
xlabel('x [mm]');
ylabel('HU');

%% Save for reuse in the following scripts
outFile = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','phantomBuilder','phantom_SOBP_slabs_2.mat');
save(outFile,'ct','cst');
