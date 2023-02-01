function stfT = OverrideStf(stf)

stfT = stf;

stfT.numOfRays = 1;

stf_ray_fields = fieldnames(stf.ray(1));
stf_ray_fields(:,2) = cell(size(stf_ray_fields));
stf_ray_fields = stf_ray_fields';
stfT.ray = struct(stf_ray_fields{:});

stfT.ray.rayPos_bev = [0 0 0];
stfT.ray.targetPoint_bev = [0 stf.SAD 0];
stfT.ray.rayPos = [0 0 0];
stfT.ray.targetPoint = stfT.ray.targetPoint_bev;
stfT.ray.energy = stf.ray(1).energy(:)';


stf_ray_rangeShifter.ID = 0;
stf_ray_rangeShifter.eqThickness = 0;
stf_ray_rangeShifter.sourceRashiDistance = 0;

stfT.ray.rangeShifter = struct();
stfT.ray.rangeShifter = repmat(stf_ray_rangeShifter,[1,size(stfT.ray.energy,2)]);


stfT.ray.focusIx = ones(1,size(stfT.ray.energy,2));
stfT.ray.numParticlesPerMU = 1000000*ones(1,size(stfT.ray.energy,2));
stfT.ray.minMU = zeros(1,size(stfT.ray.energy,2));
stfT.ray.maxMU = Inf(1,size(stfT.ray.energy,2));

stfT.numOfBixelsPerRay = size(stfT.ray.energy,2);
stfT.totalNumOfBixels = sum(stfT.numOfBixelsPerRay);
%stfT.isoCenter = [0 0 0];
end