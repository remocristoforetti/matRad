function physicalDose = matRad_loadScenariosFromMeta(saveDir, scenariosMeta, dijTemplate)

    matRad_cfg = MatRad_Config.instance();

    if ~exist('dijTemplate', 'var') || isempty(dijTemplate)

        matRad_cfg.dispInfo('No dijTemplate provided, trying to locate one ...');
        try
            load(fullfile(saveDirectory, 'dijTemplate'));
            matRad_cfg.dispInfo('done. \n');
        catch

            matRad_cfg.dispError('Unable to load dijTemplate file');
        end
    end

    % Get meta info from dij template
    nVoxels = dijTemplate.doseGrid.numOfVoxels;
    nBixels = dijTemplate.totalNumOfBixels;

    nScensToLoad = numel(scenariosMeta);
    
    
    physicalDose = arrayfun(@(scen) spalloc(nVoxels,nBixels,scen.nnzElements), scenariosMeta, 'UniformOutput',false);
    
    for scenIdx=1:nScensToLoad

        currScenarioMeta = scenariosMeta(scenIdx);
        fileName = fullfile(saveDir, currScenarioMeta.name);
        
        currDijScen = load(fileName, 'dijScenario');
        currDijScen = currDijScen.dijScenario;
        physicalDose(scenIdx) = currDijScen;
    end
end