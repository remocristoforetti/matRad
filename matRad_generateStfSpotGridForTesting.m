function stf = matRad_generateStfSpotGridForTesting(ct,cst,pln,energies,x,y)
    matRad_cfg = MatRad_Config.instance();
    
    % if ~exist('spotGridSize', 'var') || isempty(spotGridSize)
    %     spotGridSize = [1 1];
    % end

    % if ~exist('x', 'var') || isempty(x)
    %     x = zeros(spotGridSize(1)*spotGridSize(2),numel(energies));
    % end
    % 
    % if ~exist('y', 'var') || isempty(y)
    %     y = zeros(spotGridSize(1)*spotGridSize(2), numel(energies));
    % end
    % 
    % if ~exist('div', 'var') || isempty(div)
    %     div = 1;
    % end

    %[xGrid, yGrid] = meshgrid(x,y);
    
    % xGrid_V = xGrid(:);
    % yGrid_V = yGrid(:);
    xGrid_V = x;
    yGrid_V = y;
    
    nEnergies = numel(energies);
    nSpots = numel(xGrid_V);

    stf = matRad_generateStf(ct, cst,pln);
    for i=1:numel(stf)
        stf(i).ray = struct();
    % stf.ray = struct('rayPos_bev', [],...
    %     'targetPoint_bev', [],...
    %     'rayPos', [],...
    %     'targetPoint', [],...
    %     'energy', [],...
    %     'rangeShifter', [],...
    %     'focusIx', []);
        currPln = pln;
        currPln.propStf.gantryAngles = pln.propStf.gantryAngles(i);
        currPln.propStf.couchAngles  = pln.propStf.couchAngles(i);
        currPln.propStf.isoCenter    = pln.propStf.isoCenter(i,:);

        for spotIdx=1:nSpots
    
                    tmpStf = matRad_generateStfSinglePencilBeam(ct, cst,currPln,energies(1),xGrid_V(spotIdx), yGrid_V(spotIdx));
    
                    tmpRay = tmpStf.ray;
    
                    tmpFields = fieldnames(tmpRay);
                    for k=1:numel(tmpFields)
                        stf(i).ray(spotIdx).(tmpFields{k}) = tmpRay.(tmpFields{k});
                    end
                    stf(i).ray(spotIdx).energy = energies;
                    stf(i).ray(spotIdx).rangeShifter = repmat(tmpRay.rangeShifter, 1, nEnergies);
                    stf(i).ray(spotIdx).focusIx = repmat(tmpRay.focusIx, 1,nEnergies);
                    stf(i).ray(spotIdx).numParticlesPerMU = repmat(10^6, 1, nEnergies);
                    stf(i).ray(spotIdx).minMU = zeros(1,nEnergies);
                    stf(i).ray(spotIdx).maxMU = Inf(1,nEnergies);         
        end


        stf(i).numOfBixelsPerRay = nEnergies*ones(1,nSpots);
        
        stf(i).totalNumOfBixels = nEnergies*nSpots;
        stf(i).numOfRays       = nSpots;
    end

end