function ctT = OverrideCT(ct)
    ctT = ct;  
    ctT.cubeDim = [1 1 100]; %[y,x,z]
    
    ctT.cube = {ones(ctT.cubeDim)};
    
    ctT.resolution.x = 100;
    ctT.resolution.y = 100; %Resolution*nVoxels = size(phantom)
    ctT.resolution.z = 1;
    
    ctT.cubeHU = {zeros(ctT.cubeDim)};
    ctT.hlut = ct.hlut;
end