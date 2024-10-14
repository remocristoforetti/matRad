function [dij, dijTemplate] = matRad_calcParticleDoseMultipleScenarios(ct,stf,pln,cst,calcDoseDirect, saveDirectory, loadScenariosIntoDij)
% Wrapper to compute and save one scenario at a time and keepm RAM free

    matRad_cfg = MatRad_Config.instance();


    if ~exist('calcDoseDirect', 'var') || isempty(calcDoseDirect)
        calcDoseDirect = 0;
    end

    if ~exist('rmTmpFielsAfterCalc', 'var') || isempty(rmTmpFielsAfterCalc)
        rmTmpFielsAfterCalc = 0;
    end
    %saveDirectory = fullfile(matRad_cfg.matRadRoot, 'ICCRdijs', 'tmpScenarioDirectory');

    if ~exist(saveDirectory, 'dir')
        mkdir(saveDirectory);
    end

    if numel(dir(saveDirectory))>2

        matRad_cfg.dispError('temporary scenario directory contains elements, please empty the directory first');
     
    end
        
    %% compute dij scenarios
    tStart = tic;
    scenCounter = 1;
    nnzElementsV = [];
    scenIndexes = [];

    if isa(pln.multScen, 'matRad_RandomScenarios') || isa(pln.multScen, 'matRad_WorstCaseScenarios')
        totNumShiftScen     = pln.multScen.totNumShiftScen;
        numOfCtScen         = pln.multScen.numOfCtScen;
        totNumRangeScen     = pln.multScen.totNumRangeScen;
    elseif isa(pln.multScen, 'matRad_NominalScenario')
        totNumShiftScen     = 1;
        numOfCtScen         = pln.multScen.numOfCtScen;
        totNumRangeScen     = 1;
    end
    
    for shiftScenIdx=1:totNumShiftScen
        for ctScenIdx=1:numOfCtScen
            for rangeShiftScenIdx=1:totNumRangeScen

                currScenIndex = sub2ind([numOfCtScen, totNumShiftScen, totNumRangeScen], ctScenIdx, shiftScenIdx,rangeShiftScenIdx);
                %[currCTind, currShiftInd, currRangeInd] = ind2sub([pln.multScen.numOfCtScen, pln.multScen.totNumShiftScen, pln.multScen.numOfRangeShiftScen],scenCounter);
                
                if pln.multScen.scenMask(ctScenIdx, shiftScenIdx,rangeShiftScenIdx)
                    matRad_cfg.dispInfo('Computing scenario %i out of %i\n', scenCounter, pln.multScen.totNumScen);
                    currCt = ct;
                    currCt.numOfCtScen = 1;
                    currCt.cubeHU(:) = [];
                    currCt.cubeHU{1} = ct.cubeHU{ctScenIdx};
                    
                    currPln = pln;
                    currPln.multScen = matRad_multScen(currCt,'nomScen');
    
                    scenIdx = find(pln.multScen.linearMask(:,1) == ctScenIdx & ...
                                        pln.multScen.linearMask(:,2) == shiftScenIdx & ...
                                        pln.multScen.linearMask(:,3) == rangeShiftScenIdx);
                    
                    % currPln.multScen.isoShift         = pln.multScen.isoShift(scenIdx,:);
                    % currPln.multScen.relRangeShift    = pln.multScen.relRangeShift(scenIdx);
                    % currPln.multScen.absRangeShift    = pln.multScen.absRangeShift(scenIdx);
                    % currPln.multScen.maxAbsRangeShift = pln.multScen.absRangeShift(scenIdx);
                    % currPln.multScen.maxRelRangeShift = pln.multScen.relRangeShift(scenIdx);
                    % currPln.multScen.scenProb         = pln.multScen.scenProb(scenIdx);
                    currPln.multScen.isoShift         = pln.multScen.isoShift(pln.multScen.linearMask(scenIdx,2),:);
                    currPln.multScen.relRangeShift    = pln.multScen.relRangeShift(pln.multScen.linearMask(scenIdx,3));
                    currPln.multScen.absRangeShift    = pln.multScen.absRangeShift(pln.multScen.linearMask(scenIdx,3));
                    currPln.multScen.maxAbsRangeShift = pln.multScen.absRangeShift(pln.multScen.linearMask(scenIdx,3));
                    currPln.multScen.maxRelRangeShift = pln.multScen.relRangeShift(pln.multScen.linearMask(scenIdx,3));
                    currPln.multScen.scenProb         = pln.multScen.scenProb(scenIdx);
    
    
                    currCst = cst;
                    for i=1:size(cst)
                        currCst{i,4}(:) = [];
                        currCst{i,4}(1) = cst{i,4}(ctScenIdx);
                    end

                    if strcmp(pln.radiationMode, 'protons')
                        currDij = matRad_calcParticleDose(currCt,stf,currPln,currCst);
                    else
                        currDij = matRad_calcPhotonDose(currCt, stf, currPln, currCst);
                    
                    end
    
                    matRad_cfg.dispInfo('saving scenario...');
                    fileName = fullfile(saveDirectory, ['scenario_', num2str(scenIdx)]);
    
                    dijScenario = currDij.physicalDose;
                    %dijIndex = [currScenIndex];
                    % isoShift = pln.multScen.isoShift(scenIdx,:);
                    % relRangeSHift = pln.multScen.relRangeShift(scenIdx);
                    % absRangeShift = pln.multScen.absRangeShift(scenIdx);

                    isoShift = pln.multScen.isoShift(pln.multScen.linearMask(scenIdx,2),:);
                    relRangeSHift = pln.multScen.relRangeShift(pln.multScen.linearMask(scenIdx,3));
                    absRangeShift = pln.multScen.absRangeShift(pln.multScen.linearMask(scenIdx,3));

                    scenProb = pln.multScen.scenProb(scenIdx);

                    scenIndexes = [scenIndexes, currScenIndex];
                    nnzElements = nnz(currDij.physicalDose{1});
                    nnzElementsV = [nnzElementsV, nnzElements];

                    save(fileName, 'dijScenario', 'nnzElements', 'ctScenIdx', 'isoShift','relRangeSHift','absRangeShift','scenProb', '-v7.3');

                    matRad_cfg.dispInfo('done\n');

                    
                    if scenCounter ==1
                        dijTemplate = currDij;
                        dijTemplate.physicalDose = {};
                    end
                    
                    clear currDij dijScenario;

                    scenCounter = scenCounter+1;
                end
            end
        end
    end
    tEnd = toc(tStart);

    dij = dijTemplate;
    dij.doseCalcTime = tEnd;
    
    save(fullfile(saveDirectory, 'dijTemplate.mat'), 'dijTemplate');

    save(fullfile(saveDirectory, 'stf.mat'), 'stf');

    if ~exist('loadScenariosIntoDij', 'var') || isempty(loadScenariosIntoDij)
    
        loadScenariosIntoDij = false;
    end

    if loadScenariosIntoDij
        dij = matRad_loadDijScenarios(dijTemplate, pln, saveDirectory);
    end
end