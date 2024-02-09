
div = 0;

% rayPos_bev = 10*[2,0,0;...
%               1,0,0;...
%               0,0,0;...
%               0,0,1;...
%               0,0,2;...
%               0,0,3;...
%               0,0,4];
rayPos_bev = [0,0,50];

for i=1:numel(pln.propStf.gantryAngles)
    stf(i).gantryAngle = pln.propStf.gantryAngles(i);
    stf(i).couchAngle = pln.propStf.couchAngles(i);
    stf(i).bixelWidth = pln.propStf.bixelWidth;
    stf(i).radiationMode = pln.radiationMode;
    stf(i).machine = pln.machine;
    
    
    stf(i).isoCenter = pln.propStf.isoCenter(i,:);
    stf(i).numOfRays = size(rayPos_bev,1);
    machine = load('protons_generic_MCsquare.mat');
    stf(i).SAD = machine.machine.meta.SAD;
    strRS.ID = 0;
    strRS.eqThickness = 0;
    strRS.sourceRashiDistance = 0;
    rotMat = matRad_getRotationMatrix(stf(i).gantryAngle,stf(i).couchAngle);
    for j=1:stf(i).numOfRays
       stf(i).ray(j).rayPos_bev        = rayPos_bev(i,:);
       stf(i).ray(j).targetPoint_bev   = [rayPos_bev(i,1)+div*sign(rayPos_bev(i,1)),1000,rayPos_bev(i,3)+div*sign(rayPos_bev(i,1))];
       stf(i).ray(j).rayPos            = rayPos_bev(i,:)*rotMat;
       stf(i).ray(j).targetPoint       = stf(i).ray(j).targetPoint_bev*rotMat;%[1000, rayPos_bev(i,1)+div*sign(rayPos_bev(i,1)),rayPos_bev(i,3)+div*sign(rayPos_bev(i,1))];
       stf(i).ray(j).energy            = [machine.machine.data([37]).energy];
       stf(i).ray(j).focusIx           = ones(size(stf(i).ray(j).energy));
       stf(i).ray(j).numParticlesPerMU = 1000000*ones(size(stf(i).ray(j).energy));
       stf(i).ray(j).minMU             = zeros(size(stf(i).ray(j).energy));
       stf(i).ray(j).maxMU             = Inf*ones(size(stf(i).ray(j).energy));
       stf(i).ray(j).weight            = linspace(0,1,numel(stf(i).ray(j).energy));  
       stf(i).ray(j).rangeShifter      = repmat(strRS,1,numel(stf(i).ray(j).energy));
    end
    
    stf(i).sourcePoint_bev     = [0 -stf(i).SAD 0];
    
    stf(i).sourcePoint         = [0 -stf(i).SAD 0]*rotMat;
    stf(i).numOfBixelsPerRay   = numel(stf(i).ray(1).energy)*ones(size(stf(i).ray));
    
    stf(i).longitudinalSpotSpacing = 2;
    stf(i).totalNumOfBixels    = sum(stf(i).numOfBixelsPerRay);
end