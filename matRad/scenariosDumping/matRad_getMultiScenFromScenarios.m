function multScen = matRad_getMultiScenFromScenarios(saveDir,scenarioModel)
    
    matRad_cfg = MatRad_Config.instance();

    if ~exist('scenarioModel', 'var')
        scenarioModel = 'rndScen';
    end
    % Recover scenarios meta
    scenariosMeta = matRad_getMetaFromScenariosPool(saveDir);

    numDetectedCtScenarios = numel(unique([scenariosMeta.ctScenIdx]));
    numDetectedScenarios   = numel(scenariosMeta);
    %Build the multscen
    if isstruct(scenarioModel)
        if isfield(scenarioModel, 'model')
            currMultiScenMeta = scenarioModel;
        else
            matRad_cfg.dispError('Provided input scenarioModel is a struct but does not contain the field: model');
        end
    elseif ischar(scenarioModel)
        currMultiScenMeta.model = scenarioModel;
    else
        matRad_cfg.dispError('Invalid scenarioModel provided');
    end

    
    ct.numOfCtScen = numDetectedCtScenarios;

    isoShiftMat = reshape([scenariosMeta.isoShift], 3, numDetectedScenarios)';
    relShiftMat = [scenariosMeta.relRangeShift]';
    absShiftMat = [scenariosMeta.absRangeShift]';
    senProbMat  = [scenariosMeta.scenProb]';

    % Important that first field is nSamples, otherwise assigment in
    % multScen will override everything
    currMultiScenMeta.nSamples        = numDetectedScenarios./numDetectedCtScenarios;
    currMultiScenMeta.totNumShiftScen = size(unique(isoShiftMat, 'rows'),1);
    currMultiScenMeta.totNumRangeScen = size(unique([relShiftMat, absShiftMat], 'rows'),1);
    currMultiScenMeta.maxAbsRangeShift = max(absShiftMat);
    currMultiScenMeta.maxRelRangeShift = max(relShiftMat);

    currMultiScenMeta.isoShift      = isoShiftMat;
    currMultiScenMeta.relRangeShift = relShiftMat;
    currMultiScenMeta.absRangeShift = absShiftMat;
    currMultiScenMeta.scenProb      = senProbMat;

    multScen = matRad_ScenarioModel.create(currMultiScenMeta,ct);

end