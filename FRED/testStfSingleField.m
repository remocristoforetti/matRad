rayPos_bev = 20*[2,0,0;...
              1,0,0;...
              0,0,0;...
              0,0,1;...
              0,0,2;...
              0,0,3;...
              0,0,4];



stf.gantryAngle = pln.propStf.gantryAngles;
stf.couchAngle = pln.propStf.couchAngles;
stf.bixelWidth = pln.propStf.bixelWidth;
stf.radiationMode = pln.radiationMode;
stf.machine = pln.machine;

stf.SAD = 10000;
stf.isoCenter = pln.propStf.isoCenter;
stf.numOfRays = 7;
machine = load('protons_Generic.mat');
for i=1:stf.numOfRays
    stf.ray(i).rayPos_bev        = rayPos_bev(i,:);
    stf.ray(i).targetPoint_bev   = [rayPos_bev(i,1)+1*sign(rayPos_bev(i,1)),1000,rayPos_bev(i,3)+1*sign(rayPos_bev(i,1))];
    stf.ray(i).rayPos            = permute(rayPos_bev(i,:), [2,1,3]);
    stf.ray(i).targetPoint       = [1000, rayPos_bev(i,1)+1*sign(rayPos_bev(i,1)),rayPos_bev(i,3)+1*sign(rayPos_bev(i,1))];
    stf.ray(i).energy            = [machine.machine.data([30, 35, 40]).energy];
    stf.ray(i).focusIx           = ones(size(stf.ray(i).energy));
    stf.ray(i).numParticlesPerMU = 1000000*ones(size(stf.ray(i).energy));
    stf.ray(i).minMU             = zeros(size(stf.ray(i).energy));
    stf.ray(i).maxMU             = Inf*ones(size(stf.ray(i).energy));
    stf.ray(i).weight            = linspace(0,1,numel(stf.ray(i).energy));    
    stf.ray(i).rangeShifter      = [];
end

stf.sourcePoint_bev     = [0 -stf.SAD 0];

stf.sourcePoint         = [stf.SAD 0 0];
stf.numOfBixelsPerRay   = numel(stf.ray(1).energy)*ones(size(stf.ray));
stf.longitudinalSpotSpacing = 2;
stf.totalNumOfBixels    = sum(stf.numOfBixelsPerRay);