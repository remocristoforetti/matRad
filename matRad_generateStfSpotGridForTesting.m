function stf = matRad_generateStfSpotGridForTesting(ct,cst,pln,energies,x,y)
    matRad_cfg = MatRad_Config.instance();
    
    % if ~exist('spotGridSize', 'var') || isempty(spotGridSize)
    %     spotGridSize = [1 1];
    % end

    if ~exist('x', 'var') || isempty(x)
        x = zeros(spotGridSize(1)*spotGridSize(2),numel(energies));
    end
    
    if ~exist('y', 'var') || isempty(y)
        y = zeros(spotGridSize(1)*spotGridSize(2), numel(energies));
    end
    
    if ~exist('div', 'var') || isempty(div)
        div = 1;
    end

    [xGrid, yGrid] = meshgrid(x,y);

    xGrid_V = xGrid(:);
    yGrid_V = yGrid(:);

    nEnergies = numel(energies);
    nSpots = numel(xGrid_V);

    stf = matRad_generateStf(ct, cst,pln);
    stf.ray = struct();
    % stf.ray = struct('rayPos_bev', [],...
    %     'targetPoint_bev', [],...
    %     'rayPos', [],...
    %     'targetPoint', [],...
    %     'energy', [],...
    %     'rangeShifter', [],...
    %     'focusIx', []);
    
    for spotIdx=1:nSpots

                tmpStf = matRad_generateStfSinglePencilBeam(ct, cst,pln,energies(1),xGrid_V(spotIdx), yGrid_V(spotIdx));

                tmpRay = tmpStf.ray;

                tmpFields = fieldnames(tmpRay);
                for k=1:numel(tmpFields)
                    stf.ray(spotIdx).(tmpFields{k}) = tmpRay.(tmpFields{k});
                end
                stf.ray(spotIdx).energy = energies;
                stf.ray(spotIdx).rangeShifter = repmat(tmpRay.rangeShifter, 1, nEnergies);
                stf.ray(spotIdx).focusIx = repmat(tmpRay.focusIx, 1,nEnergies);
                stf.ray(spotIdx).numParticlesPerMU = repmat(10^6, 1, nEnergies);
                stf.ray(spotIdx).minMU = zeros(1,nEnergies);
                stf.ray(spotIdx).maxMU = Inf(1,nEnergies);         
    end


    stf.numOfBixelsPerRay = nEnergies*ones(1,nSpots);
    
    stf.totalNumOfBixels = nEnergies*nSpots;
    stf.numOfRays       = nSpots;
end