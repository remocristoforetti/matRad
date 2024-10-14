function [dij,currMultiscen,loadedScenarios] = matRad_loadDijScenarios(ct,saveDirectory, loadScenMode, nScensToLoad, verboseLevel, scenarioCombination, wcSigma)
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
    

    if ~exist('verboseLevel', 'var') || isempty(verboseLevel)
        verboseLevel = 1;
    end

    if ~exist('loadScenMode', 'var') || isempty(loadScenMode)
         loadScenMode = 'all';
    end

    if ~exist('scenarioCombination', 'var') || isempty(scenarioCombination)
         scenarioCombination = 'rnd';
    end

    if ~exist('wcSigma', 'var') || isempty(wcSigma)
        wcSigma = 1;
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
            case 'none'
                nAllScens = sum(arrayfun(@(file) contains(file.name, 'scenario'), dirFiles(3:end)));%sum(arrayfun(@(file) strcmp(file.name(1:8), 'scenario'), dirFiles(3:end)));
                scensToLoadMeta = [1:nAllScens];
                loadScenarios = false;
        end
    else
        if isnumeric(loadScenMode)
            nAllScens = numel(loadScenMode);
            if ~isrow(loadScenMode)
                loadScenMode = loadScenMode';
            end
            scensToLoadMeta = loadScenMode;

        end
    end

    if exist('nScensToLoad', 'var') && ischar(nScensToLoad) && strcmp(nScensToLoad, 'none')
        loadScenarios = false;
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
    
    % This is the original multiscenario combination (random or simple worst case for now)
    switch scenarioCombination
        case 'rnd'
            % Valid for random generated scenarios
            %This is number of error scenarios only (total number of scenarios/nCtScenario)
            nErrorScenarios = totNumOfShiftScen;
            if ~exist('nScensToLoad', 'var') || isempty(nScensToLoad)
                %load all the error scenarios
                errorsToLoad = [1:nErrorScenarios];
            else
                % 
                % if nScensToLoad <= nErrorScenarios
                % 
                %     %Peak randomly which scenarios to keep
                %     errorsToLoad = sort(randperm(nErrorScenarios,nScensToLoad));
                % else
                %     errorsToLoad = [1:nErrorScenarios];
                % end
                if ischar(nScensToLoad)
                    errorsToLoad = [1:nErrorScenarios];
                else
                    errorsToLoad = nScensToLoad;
                end
            end
        
            if nAllScens>1
                if totNumOfRangeShift>=1
                    currMultiscen = matRad_multScen(ct, 'rndScen');
                    currMultiscen.nSamples = numel(errorsToLoad);
                    currMultiscen.isoShift = shiftScenMat(errorsToLoad,:);%unique(shiftScenMat, 'rows', 'stable');
                    currMultiscen.relRangeShift = relAbsErrorMat(errorsToLoad,1);
                    currMultiscen.absRangeShift = relAbsErrorMat(errorsToLoad,2);
                    currMultiscen.maxAbsRangeShift = max(relAbsErrorMat(errorsToLoad,2));
                    currMultiscen.maxRelRangeShift = max(relAbsErrorMat(errorsToLoad,1));
                else            
                    currMultiscen = matRad_multScen(ct, 'rndScen_shiftOnly');          
                    currMultiscen.nSamples = numel(errorsToLoad);
                    currMultiscen.isoShift = shiftScenMat(errorsToLoad,:);%unique(shiftScenMat, 'rows', 'stable');
                end
            else
                    currMultiscen = matRad_multScen(ct, 'rndScen');
                    currMultiscen.nSamples = numel(errorsToLoad);
                    currMultiscen.isoShift = shiftScenMat(errorsToLoad,:);%unique(shiftScenMat, 'rows', 'stable');
                    currMultiscen.relRangeShift = relAbsErrorMat(errorsToLoad,1);
                    currMultiscen.absRangeShift = relAbsErrorMat(errorsToLoad,2);
                    currMultiscen.maxAbsRangeShift = max(relAbsErrorMat(errorsToLoad,2));
                    currMultiscen.maxRelRangeShift = max(relAbsErrorMat(errorsToLoad,1));
            end
            
            % Keep only the scenarios that have errors equal to the randomly sampled ones 
            loadScens = ismember(shiftScenAll, currMultiscen.isoShift,'rows') & ismember(relAbsScenAll, [currMultiscen.relRangeShift, currMultiscen.absRangeShift], 'rows');
            
            % for random scenarios only
        
        case 'wc'

            nErrorScenarios = totNumOfShiftScen + totNumOfRangeShift -1; % -1 because nominal scenario is included
            
         
            %load all the error scenarios
            errorsToLoad = [1:nErrorScenarios];
     
            % if numel(errorsToLoad) ~= nErrorScenarios
            %     matRad_cfg.dispError('For wc scenarios need to pick all scenarios');
            % end

            currMultiscen = matRad_multScen(ct, 'wcScen');
            currMultiscen.wcSigma = wcSigma;
            %currMultiscen.relRangeShift = [scenarios.relRangeSHift]';
            %currMultiscen.absRangeShift = [scenarios.absRangeShift]';
            %currMultiscen.maxAbsRangeShift = max([scenarios.absRangeShift]);
            %currMultiscen.maxRelRangeShift = max([scenarios.relRangeSHift]);

            %for i=1:numel(scenarios)
            %    currMultiscen.isoShift(i,:) = scenarios(i).isoShift;
            %end

            loadScens = ones(numel(scenarios),1);

        case 'wc_shift'
             nErrorScenarios = totNumOfShiftScen + totNumOfRangeShift -1; % -1 because nominal scenario is included
            
         
            %load all the error scenarios
            errorsToLoad = [1:nErrorScenarios];


            currMultiscen = matRad_multScen(ct, 'wcScen');
            currMultiscen.wcSigma = wcSigma;
            currMultiscen.combinations = 'shift';
            currMultiscen.combineRange = true;

            loadScens = ones(numel(scenarios),1);
    end

    
     numOfShiftScen = currMultiscen.totNumShiftScen;
     numOfRangeShift = currMultiscen.totNumRangeScen;
  

    %dij.physicalDose = cell(numOfCtScen, numOfShiftScen, numOfRangeShift);
    dij.physicalDose = cell(currMultiscen.numOfCtScen, numOfShiftScen, numOfRangeShift);
    if isnumeric(loadScenMode) && currMultiscen.totNumScen ~= numel(loadScenMode)
        matRad_cfg.dispWarning('Loading scenarios: but inconsistent multipleScenario output provided');
    end

     %if loadScenarios    
        for scenIdx=1:nAllScens    
            if loadScens(scenIdx)
                currCTidx = scenarios(scenIdx).ctScenIdx;

                switch scenarioCombination
                    case 'rnd'
                        currShiftScenIdx = find(ismember(currMultiscen.isoShift,scenarios(scenIdx).isoShift, 'rows'));
                
                        if numOfRangeShift>1
                            currAbsRelShiftdx = find(ismember([currMultiscen.relRangeShift, currMultiscen.absRangeShift], [scenarios(scenIdx).relRangeSHift,scenarios(scenIdx).absRangeShift], 'rows'));
                        else
                            currAbsRelShiftdx = 1;
            
                        end
                        
                        scenarios(scenIdx).scenarioIndexDij = sub2ind([currMultiscen.numOfCtScen, numOfShiftScen, numOfRangeShift], currCTidx, currShiftScenIdx, currAbsRelShiftdx);
                   
                        %currMultiscen.scenMask(scenarios(scenIdx).scenarioIndexDij) = true;
                        linIndex = find(ismember(currMultiscen.linearMask, [currCTidx,currShiftScenIdx,currAbsRelShiftdx], 'rows'));
                        currMultiscen.scenProb(linIndex) = scenarios(scenIdx).scenProb;
                        %dij.physicalDose{scenarios(scenIdx).scenarioIndexDij} = spalloc(dij.doseGrid.numOfVoxels, dij.totalNumOfBixels, scenarios(scenIdx).nnzElements);

                    case {'wc', 'wc_shift'}

                        currShiftScen = scenarios(scenIdx).isoShift;
                        currRelRange  = scenarios(scenIdx).relRangeSHift;
                        currAbsRange  = scenarios(scenIdx).absRangeShift;

                        if currAbsRange == 0 && currRelRange == 0
                            % these are shift scenarios
                            currRangeErrorLinearIdx = 1;
                            if all(currShiftScen == [0,0,0])
                                % this is teh nominal scenario
                                shiftScenLinearIdx = 1;
                            else
                                shiftScenLinearIdx = find(ismember(currMultiscen.isoShift, currShiftScen, 'rows'));
                            end
                        else
                            matRad_cfg.dispWarning('!!! This setup only works for 1 single CT scenario');
                            % These are range error scenarios
                            shiftScenLinearIdx = 1;
                            rangeErrorOnlyScenarios = [currMultiscen.relRangeShift,currMultiscen.absRangeShift];
                            rangeErrorOnlyScenarios(find(ismember(rangeErrorOnlyScenarios, [0,0], 'rows')),:) = [];
                            currRangeErrorLinearIdx = find(ismember(rangeErrorOnlyScenarios, [currRelRange,currAbsRange], 'rows'));
                        end

                        scenLinearIdx = find(ismember(currMultiscen.linearMask, [currCTidx, shiftScenLinearIdx, currRangeErrorLinearIdx], 'rows'));
                      
                        scenarios(scenIdx).scenarioIndexDij = sub2ind([currMultiscen.numOfCtScen, numOfShiftScen, numOfRangeShift], currMultiscen.linearMask(scenLinearIdx,1), currMultiscen.linearMask(scenLinearIdx,2), currMultiscen.linearMask(scenLinearIdx,3));

                       
                        linIndex = scenIdx;
                        currMultiscen.scenProb(linIndex) = scenarios(scenIdx).scenProb;
                end
            end
        end

        % % This needs to be changed later
        % (18/9/2024)
        % if nAllScens == 1
        %     scenarios(1).scenarioIndexDij = 1;
        % end
        
        if loadScenarios
            for scenIdx=1:nAllScens
                dij.physicalDose{scenarios(scenIdx).scenarioIndexDij} = spalloc(dij.doseGrid.numOfVoxels, dij.totalNumOfBixels, scenarios(scenIdx).nnzElements);
            end
        end

        % Load the scenarios
        stringLenght = 0;
        scenCounter = 1;
        if loadScenarios
            for scenIdx=1:nAllScens
                if loadScens(scenIdx)
                    if verboseLevel
                        fprintf(repmat('\b',1,stringLenght));
                        stringLenght = fprintf('Loading scenarios: %u/%u\n', scenCounter, sum(loadScens));
                    end
                    fileIndex = find(strcmp({dirFiles.name},scenarios(scenIdx).name));
            
                    currDijScen = load(fullfile(dirFiles(fileIndex).folder, dirFiles(fileIndex).name), 'dijScenario');
                    dij.physicalDose(scenarios(scenIdx).scenarioIndexDij) = currDijScen.dijScenario;
                    scenCounter = scenCounter+1;
                end
         
            end
        end
    loadedScenarios = {scenarios(loadScens).name};
end