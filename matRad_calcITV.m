function [cstITV, ctITV] = matRad_calcITV(cst, ct)

    %[cstITV, ctITV] = matRad_calcITV(cst, ct)
    
    matRad_cfg = MatRad_Config.instance();


    nStructs = size(cst,1);
    nCtScenarios = ct.numOfCtScen;
    targetIdx = find(strcmp(cst(:,3), 'TARGET'));

    if isempty(targetIdx)
        matRad_cfg.dispError('No TARGET structure found in cst');
    elseif numel(targetIdx)>1
        matRad_cfg.dispWarning('Multiple TARGET structures detected, selecting: %s as reference', cst{targetIdx(1),3});
        targetIdx = targetIdx(1);
    end

    
    % collect all voxels in ct scenarios
    allVoxels = arrayfun(@(scenStruct) scenStruct{1}', cst{targetIdx,4}, 'UniformOutput',false);
    voxelIdxITV = unique([allVoxels{:}])';

    itvCstIdx = nStructs+1;
    cstITV = cst;
    cstITV{itvCstIdx,1} = itvCstIdx-1;
    cstITV{itvCstIdx,2} = 'ITV';
    cstITV{itvCstIdx,3} = 'TARGET';
    cstITV{itvCstIdx,4} = repmat({voxelIdxITV}, 1, nCtScenarios);
    cstITV{itvCstIdx,5} = struct('Priority', 1, ...
                              'alphaX', 0.1, ...
                              'betaX', 0.05, ...
                              'Visible', 1, ...
                              'visibleColor', [1 1 0]);
    cstITV{itvCstIdx,6} = [];

    % override the ct values
    ctITV = ct;
    % for ctIdx=1:ct.numOfCtScen
    %     ctITV.cubeHU{ctIdx}(voxelIdxITV) = 0;
    % end

    ctITV.numOfCtScen = 1;
    ctITV.cubeHU = ctITV.cubeHU(1);
    ctITV.cubeHU{1}(voxelIdxITV) = 0;
    

end