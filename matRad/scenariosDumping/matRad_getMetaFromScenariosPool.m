function [scenariosMeta, dijTemplate] =  matRad_getMetaFromScenariosPool(saveDir,scenIdxToInclude)
    
    matRad_cfg = MatRad_Config.instance();

    if ~exist(saveDir, 'dir')
        matRad_cfg.dispError('Could not find directory: %s', saveDir);
    end

    try
        load(fullfile(saveDir, 'dijTemplate.mat'),'dijTemplate');
        matRad_cfg.dispInfo('done. \n');
    catch
        matRad_cfg.dispWarning('Unable to load dijTemplate file');
    end

    % Detect files
    dirFiles = dir(saveDir);

    numDetectedScenarios = sum(arrayfun(@(file) contains(file.name, 'scenario'), dirFiles(3:end)));

    if ~exist('scenIdxToInclude', 'var') || isempty(scenIdxToInclude)
        scenIdxToInclude = [1:numDetectedScenarios];
    end

    % Load the meta information for all the detected scenarios
    scenariosMeta = [];
    for scenIdx=scenIdxToInclude

       fileIndex = find(strcmp({dirFiles.name},['scenario_', num2str(scenIdx), '.mat']));
       if ~isempty(fileIndex)
           fileVarList = whos('-file', fullfile(dirFiles(fileIndex).folder,dirFiles(fileIndex).name));
        
           includeVariables = ~strcmp({fileVarList.name}, 'dijScenario');
        
           tmpStruct = struct();
           tmpStruct = load(fullfile(dirFiles(fileIndex).folder, dirFiles(fileIndex).name), fileVarList(includeVariables).name);
           tmpStruct.name = dirFiles(fileIndex).name;
           scenariosMeta = [scenariosMeta,tmpStruct];
       else
            continue;
       end
        
   end
end