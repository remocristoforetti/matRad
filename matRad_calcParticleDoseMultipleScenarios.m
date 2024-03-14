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
    for shiftScenIdx=1:pln.multScen.totNumShiftScen
        for ctScenIdx=1:pln.multScen.numOfCtScen
            for rangeShiftScenIdx=1:pln.multScen.totNumRangeScen

                currScenIndex = sub2ind([pln.multScen.numOfCtScen, pln.multScen.totNumShiftScen, pln.multScen.numOfRangeShiftScen], ctScenIdx, shiftScenIdx,rangeShiftScenIdx);
                %[currCTind, currShiftInd, currRangeInd] = ind2sub([pln.multScen.numOfCtScen, pln.multScen.totNumShiftScen, pln.multScen.numOfRangeShiftScen],scenCounter);
                
                if pln.multScen.scenMask(ctScenIdx, shiftScenIdx,rangeShiftScenIdx)
                    matRad_cfg.dispInfo('Computing scenario %i out of %i\n', scenCounter, pln.multScen.totNumScen);
                    currCt = ct;
                    currCt.numOfCtScen = 1;
                    currCt.cubeHU(:) = [];
                    currCt.cubeHU{1} = ct.cubeHU{ctScenIdx};
                    
                    currPln = pln;
                    currPln.multScen = matRad_multScen(currCt,'nomScen');
    
                    currPln.multScen.isoShift = pln.multScen.isoShift(shiftScenIdx,:);
                    currPln.multScen.relRangeShift = pln.multScen.relRangeShift(rangeShiftScenIdx);
                    currPln.multScen.absRangeShift = pln.multScen.absRangeShift(rangeShiftScenIdx);
                    currPln.multScen.maxAbsRangeShift = pln.multScen.absRangeShift(rangeShiftScenIdx);
                    currPln.multScen.maxRelRangeShift = pln.multScen.relRangeShift(rangeShiftScenIdx);
                    currPln.multScen.scenProb         = pln.multScen.scenProb(scenCounter);
    
    
                    currCst = cst;
                    for i=1:size(cst)
                        currCst{i,4}(:) = [];
                        currCst{i,4}(1) = cst{i,4}(ctScenIdx);
                    end
                    currDij = matRad_calcParticleDose(currCt,stf,currPln,currCst);
    
    
                    matRad_cfg.dispInfo('saving scenario...');
                    fileName = fullfile(saveDirectory, ['scenario_', num2str(scenCounter)]);
    
                    dijScenario = currDij.physicalDose;
                    %dijIndex = [currScenIndex];
                    isoShift = pln.multScen.isoShift(shiftScenIdx,:);
                    relRangeSHift = pln.multScen.relRangeShift(rangeShiftScenIdx);
                    absRangeShift = pln.multScen.absRangeShift(rangeShiftScenIdx);
                    scenProb = pln.multScen.scenProb(scenCounter);

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
    if loadScenariosIntoDij
        dij = matRad_loadDijScenarios(dijTemplate, pln, saveDirectory);
    end
end