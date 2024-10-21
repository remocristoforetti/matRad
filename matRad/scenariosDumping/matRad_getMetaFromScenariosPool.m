function scenariosMeta =  matRad_getMetaFromScenariosPool(saveDir)
    
    % Detect files
    dirFiles = dir(saveDir);

    numDetectedScenarios = sum(arrayfun(@(file) contains(file.name, 'scenario'), dirFiles(3:end)));

    % Load the meta information for all the detected scenarios
    scenariosMeta = [];
    for scenIdx=1:numDetectedScenarios

       fileIndex = find(strcmp({dirFiles.name},['scenario_', num2str(scenIdx), '.mat']));
       fileVarList = whos('-file', fullfile(dirFiles(fileIndex).folder,dirFiles(fileIndex).name));
        
       includeVariables = ~strcmp({fileVarList.name}, 'dijScenario');
        
       tmpStruct = struct();
       tmpStruct = load(fullfile(dirFiles(fileIndex).folder, dirFiles(fileIndex).name), fileVarList(includeVariables).name);
       tmpStruct.name = dirFiles(fileIndex).name;
       scenariosMeta = [scenariosMeta,tmpStruct];
        
   end
end