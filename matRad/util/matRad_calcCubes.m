function resultGUI = matRad_calcCubes(w,dij,scenNum,boolInterpolate, quantities)

    matRad_cfg = MatRad_Config.instance();

    if  ~exist('scenNum', 'var') || isempty(scenNum)
        scenNum = 1;
    end
    
    if ~exist('boolInterpolate', 'var') || isempty(boolInterpolate)
        boolInterpolate = true;
    end
    
    resultGUI.w = w;
    
    if isfield(dij,'numParticlesPerMU')
        resultGUI.MU = (w.*1e6) ./ dij.numParticlesPerMU;
    end
    
    % get bixel - beam correspondence  
    for i = 1:dij.numOfBeams
        beamInfo(i).suffix = ['_beam', num2str(i)];
        beamInfo(i).logIx  = (dij.beamNum == i);
    end
    beamInfo(dij.numOfBeams+1).suffix = '';
    beamInfo(dij.numOfBeams+1).logIx  = true(size(resultGUI.w,1),1);
    
    if isfield(dij, 'physicalDose') && ~isempty(dij.physicalDose)
        [ctScen,~] = ind2sub(size(dij.physicalDose),scenNum);
    elseif isfield(dij, 'physicalDoseExp') && ~isempty(dij.physicalDoseExp)
        [ctScen,~] = ind2sub(size(dij.physicalDoseExp),scenNum);
    else
        ctScen = 1;
    end

    defaultQuantities = {'physicalDose', 'LET', 'ConstRBExDose', 'RBExDose'};

    if ~exist('quantities', 'var') || isempty(quantities)
        quantities = defaultQuantities;
    else
        quantities = [quantities, defaultQuantities];
    end

    %% Compute
    for quantityIdx=1:numel(quantities)
        try
            currQuantity = quantities{quantityIdx};
            qtInstance = matRad_BackProjection.getQuantityInstanceFromName(currQuantity);
            qtInstance.useScenarios = scenNum;
            
            if isa(qtInstance, 'matRad_DistributionQuantity')
                resultGUI.(currQuantity)  = reshape(qtInstance.computeQuantity(dij,scenNum, w),dij.doseGrid.dimensions);
            end
        catch

        end
    end

    %% Interpolate
    if boolInterpolate
        if isfield(dij,'ctGrid') && any(dij.ctGrid.dimensions~=dij.doseGrid.dimensions)
            myFields = fieldnames(resultGUI);
            for i = 1:numel(myFields)
                if numel(resultGUI.(myFields{i})) == dij.doseGrid.numOfVoxels
        
                    % interpolate!
                    resultGUI.(myFields{i}) = matRad_interp3(dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z, ...
                        resultGUI.(myFields{i}), ...
                        dij.ctGrid.x,dij.ctGrid.y',dij.ctGrid.z,'linear',0);
        
                end
            end
        end
    end
end