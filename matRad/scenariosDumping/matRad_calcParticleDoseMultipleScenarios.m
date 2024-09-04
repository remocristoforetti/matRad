function [dij, dijTemplate] = matRad_calcParticleDoseMultipleScenarios(ct,stf,pln,cst,calcDoseDirect, saveDirectory)
% Wrapper to compute and save one scenario at a time and keepm RAM free

    matRad_cfg = MatRad_Config.instance();


    if ~exist('calcDoseDirect', 'var') || isempty(calcDoseDirect)
        calcDoseDirect = 0;
    end

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
    
    % Maximum effort
    for ctScenIdx=1:numOfCtScen
        for shiftScenIdx=1:totNumShiftScen    
            for rangeShiftScenIdx=1:totNumRangeScen

                %currScenIndex = sub2ind([numOfCtScen, totNumShiftScen, totNumRangeScen], ctScenIdx, shiftScenIdx,rangeShiftScenIdx);
                [~,currScenIndex] = intersect(pln.multScen.linearMask, [ctScenIdx,shiftScenIdx, rangeShiftScenIdx], 'rows');
                
                if pln.multScen.scenMask(ctScenIdx, shiftScenIdx,rangeShiftScenIdx)
                    
                    matRad_cfg.dispInfo('Computing scenario %i out of %i\n', scenCounter, pln.multScen.totNumScen);
                    currCt = ct;
                    currCt.numOfCtScen = 1;
                    currCt.cubeHU(:) = [];
                    currCt.cubeHU{1} = ct.cubeHU{ctScenIdx};
                    
                    currPln = pln;


                    currMultiScenMeta.model = 'nomScen';
                    currMultiScenMeta.isoShift      = pln.multScen.isoShift(currScenIndex,:);
                    currMultiScenMeta.relRangeShift = pln.multScen.relRangeShift(currScenIndex);
                    currMultiScenMeta.absRangeShift = pln.multScen.absRangeShift(currScenIndex);
                    currMultiScenMeta.scenProb      = pln.multScen.scenProb(currScenIndex);

                    currPln.multScen = matRad_ScenarioModel.create(currMultiScenMeta,currCt);
    
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
                    fileName = fullfile(saveDirectory, ['scenario_', num2str(currScenIndex)]);
    
                    dijScenario = currDij.physicalDose;
                   
                    isoShift = currMultiScenMeta.isoShift;
                    relRangeShift = currMultiScenMeta.relRangeShift;
                    absRangeShift = currMultiScenMeta.absRangeShift;

                    scenProb = currMultiScenMeta.scenProb;

                    scenIndexes = [scenIndexes, currScenIndex];
                    nnzElements = nnz(currDij.physicalDose{1});
                    nnzElementsV = [nnzElementsV, nnzElements];

                    save(fileName, 'dijScenario', 'nnzElements', 'ctScenIdx', 'isoShift','relRangeShift','absRangeShift','scenProb', '-v7.3');

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

end