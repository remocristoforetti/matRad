function [dij,currMultiscen,loadedScenarios] = matRad_loadDijScenarios(ct,saveDirectory, loadScenMode, nScensToLoad)
%    matRad_loadDijScenarios(ct,saveDirectloadory, loadScenMode, nScensToLoad)
%
%   loasScenMode = all, rnd, [1 2 3]   
%
%
%
    matRad_cfg = MatRad_Config.instance();
    %This is just boolean to set to zero only in case function is used to
    %recover the multiScen class for example
    
    loadScenarios = true;

    try
    
        load(fullfile(saveDirectory, 'dijTemplate'));
    catch
        matRad_cfg.dispError('Unable to load dijTemplate file');
    end

    dij = dijTemplate;
    %Get filenames in the directory

    dirFiles = dir(saveDirectory);
        
    if ~exist('loadScenMode', 'var') || isempty(loadScenMode)
         loadScenMode = 'all';
    end
    %nAllScens = sum(arrayfun(@(file) strcmp(file.name(1:8), 'scenario'), dirFiles(3:end)));
    if ischar(loadScenMode)
        switch loadScenMode
            case 'all'
                nAllScens = sum(arrayfun(@(file) contains(file.name, 'scenario'), dirFiles(3:end)));%sum(arrayfun(@(file) strcmp(file.name(1:8), 'scenario'), dirFiles(3:end)));
                scensToLoadMeta = [1:nAllScens];
            case 'rnd'
                if ~exist('nScensToLoad', 'var') || isempty(nScensToLoad)
                    matRad_cfg.dispError('rnd scneario loading, specify numebr of error scenarios to be sampled');
                end
                nAllScens = sum(arrayfun(@(file) contains(file.name, 'scenario'), dirFiles(3:end)));%sum(arrayfun(@(file) strcmp(file.name(1:8), 'scenario'), dirFiles(3:end)));
                scensToLoadMeta = [1:nAllScens];
        end
    else
        if isnumeric(loadScenMode)
            nAllScens = numel(loadScenMode);
            if ~isrow(loadScenMode)
                loadScenMode = loadScenMode';
            end
            scensToLoadMeta = loadScenMode;

            if exist('nScensToLoad', 'var') && ischar(nScensToLoad) && strcmp(nScensToLoad, 'none')
                loadScenarios = false;
            end
        end
    end


    scenarios = [];
    for scenIdx=scensToLoadMeta

        fileIndex = find(strcmp({dirFiles.name},['scenario_', num2str(scenIdx), '.mat']));

        fileVarList = whos('-file', fullfile(dirFiles(fileIndex).folder,dirFiles(fileIndex).name));
        
        includeVariables = ~strcmp({fileVarList.name}, 'dijScenario');
        


        tmpStruct = struct();
        tmpStruct = load(fullfile(dirFiles(fileIndex).folder, dirFiles(fileIndex).name), fileVarList(includeVariables).name);
        tmpStruct.name = dirFiles(fileIndex).name;
        scenarios = [scenarios,tmpStruct];%load(dirFiles(fileIndex).name, fileVarList(includeVariables).name);
        
    end


    %Allocate Space
    numOfCtScen = numel(unique([scenarios.ctScenIdx]));

    shiftScenAll = reshape([scenarios.isoShift], 3, [])';
    
    shiftScenMat = unique(shiftScenAll, 'rows', 'stable');
    totNumOfShiftScen = size(shiftScenMat,1);

    relAbsScenAll = [[scenarios.relRangeSHift]', [scenarios.absRangeShift]'];
    relAbsErrorMat = unique(relAbsScenAll, 'rows', 'stable');
    totNumOfRangeShift = size(relAbsErrorMat,1);
    
    % Valid for random generated scenarios

    %This is number of error scenarios only (total number of scenarios/nCtScenario)
    nErrorScenarios = totNumOfShiftScen;
    if ~exist('nScensToLoad', 'var') || isempty(nScensToLoad)
        %load all the error scenarios
        errorsToLoad = [1:nErrorScenarios];
    else

        if nScensToLoad <= nErrorScenarios
 
            %Peak randomly which scenarios to keep
            errorsToLoad = sort(randperm(nErrorScenarios,nScensToLoad));
        else
            errorsToLoad = [1:nErrorScenarios];
        end

    end

    currMultiscen = matRad_multScen(ct, 'rndScen');
    currMultiscen.nSamples = numel(errorsToLoad);
    currMultiscen.isoShift = shiftScenMat(errorsToLoad,:);%unique(shiftScenMat, 'rows', 'stable');
    currMultiscen.relRangeShift = relAbsErrorMat(errorsToLoad,1);
    currMultiscen.absRangeShift = relAbsErrorMat(errorsToLoad,2);
    currMultiscen.maxAbsRangeShift = max(relAbsErrorMat(errorsToLoad,2));
    currMultiscen.maxRelRangeShift = max(relAbsErrorMat(errorsToLoad,1));

    % Keep only the scenarios that have errors equal to the randomly sampled ones 
    loadScens = ismember(shiftScenAll, currMultiscen.isoShift,'rows') & ismember(relAbsScenAll, [currMultiscen.relRangeShift, currMultiscen.absRangeShift], 'rows');
    
    % for random scenarios only
    numOfShiftScen = numel(errorsToLoad);
    numOfRangeShift = numel(errorsToLoad);

    dij.physicalDose = cell(numOfCtScen, numOfShiftScen, numOfRangeShift);

    if isnumeric(loadScenMode) && currMultiscen.totNumScen ~= numel(loadScenMode)
        matRad_cfg.dispWarning('Loading scenarios: %d but inconsistent multipleScenario output provided', loadScenMode);
    end
    
        %nScens = currMultiscen.totNumScen;
    for scenIdx=1:nAllScens

        if loadScens(scenIdx)
            currCTidx = scenarios(scenIdx).ctScenIdx;
            currShiftScenIdx = find(ismember(currMultiscen.isoShift,scenarios(scenIdx).isoShift, 'rows'));
            currAbsRelShiftdx = find(ismember([currMultiscen.relRangeShift, currMultiscen.absRangeShift], [scenarios(scenIdx).relRangeSHift,scenarios(scenIdx).absRangeShift], 'rows'));
            
            scenarios(scenIdx).scenarioIndexDij = sub2ind([numOfCtScen, numOfShiftScen, numOfRangeShift], currCTidx, currShiftScenIdx, currAbsRelShiftdx);
       
            %currMultiscen.scenMask(scenarios(scenIdx).scenarioIndexDij) = true;
            linIndex = find(ismember(currMultiscen.linearMask, [currCTidx,currShiftScenIdx,currAbsRelShiftdx], 'rows'));
            currMultiscen.scenProb(linIndex) = scenarios(scenIdx).scenProb;
            dij.physicalDose{scenarios(scenIdx).scenarioIndexDij} = spalloc(dij.doseGrid.numOfVoxels, dij.totalNumOfBixels, scenarios(scenIdx).nnzElements);
        end

    end

    % Load the scenarios
    stringLenght = 0;
    scenCounter = 1;
    if loadScenarios
        for scenIdx=1:nAllScens
            if loadScens(scenIdx)
                fprintf(repmat('\b',1,stringLenght));
                stringLenght = fprintf('Loading scenarios: %u/%u\n', scenCounter, sum(loadScens));
        
                fileIndex = find(strcmp({dirFiles.name},scenarios(scenIdx).name));
        
                currDijScen = load(fullfile(dirFiles(fileIndex).folder, dirFiles(fileIndex).name), 'dijScenario');
                dij.physicalDose(scenarios(scenIdx).scenarioIndexDij) = currDijScen.dijScenario;
                scenCounter = scenCounter+1;
            end
     
        end
    end

    loadedScenarios = {scenarios(loadScens).name};
end