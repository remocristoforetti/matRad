%% Script 1: Build an elongated box phantom for single-PB / small-SOBP demos
%
% Uses matRad_PhantomBuilder to construct a homogeneous water phantom,
% elongated along x, with a small cubic target offset from the phantom
% center. The elongation along x provides enough water "run-up" for a
% particle beam entering from the -x face, and a safety margin behind the
% target so the beam can come to rest inside the phantom.
%
% Later scripts in this series reuse this ct/cst to shoot a single pencil
% beam and a small spread-out Bragg peak (SOBP) with particle therapy.

%% set matRad runtime configuration
clear;
matRad_rc; %If this throws an error, run it from the matRad root directory first to set the paths
matRad_cfg = MatRad_Config.instance();

%% Phantom / target geometry
ctDim        = [200 50 50];    % [x y z] number of voxels
ctResolution = [1 1 1];        % [x y z] mm per voxel -> phantom spans 200 x 50 x 50 mm^3
numOfCtScen  = 1;              % single (static) ct scenario, no motion

phantomHU    = 0;

targetDim    = [30 30 30];     % [x y z] mm, small cubic target
targetOffset = [60 0 0];       % [x y z] mm, shifted from the phantom center along x
                                % (beam will later enter from the -x face and travel in +
targetHU     = phantomHU;

insertDim    = [30 25 25];
insertOffset = [0 15 0];
insertHU     = phantomHU + 200;

builder = matRad_PhantomBuilder(ctDim,ctResolution,numOfCtScen);

%% Define the dose objective for the target
% matRad_SquaredDeviation(penalty,dRef) penalizes deviations from dRef.
% dRef is given as *total* physical dose in Gy; matRad_fluenceOptimization
% divides it by pln.numOfFractions internally, so this value should match
% the total prescribed dose used by the optimization scripts that follow.
targetObjective = struct(DoseObjectives.matRad_SquaredDeviation(1,60));

%% Add the volumes of interest
% The whole cube is first filled with water background (HU = 0) via a box
% OAR spanning the full ct grid, then the target box is carved out on top.
% Note: for overlapping VOIs, matRad_PhantomBuilder.getctcst() assigns HU
% in reverse definition order, so the *first-added* VOI's HU wins the
% overlap; the same VOI also gets the higher overlap *priority* used
% during optimization (cst{:,5}.Priority = row index, lower = higher
% priority). Here it makes no practical difference, since target and
% background share the same HU (water), and the background carries no
% objective, so it can never claim the target's voxels during
% optimization (see matRad_setOverlapPriorities).
builder.addBoxTarget('TARGET',targetDim,'HU',targetHU,'offset',targetOffset,'objectives',targetObjective);
builder.addBoxOAR('Insert',insertDim,'HU',insertHU, 'offset', insertOffset);
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
outFile = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','phantomBuilder','phantom_insert.mat');
save(outFile,'ct','cst');
