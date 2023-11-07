voiIx = [];
for i = 1:size(cst,1)
    if ~isempty(cst{i,6})
        voiIx = [voiIx i];
    end
end

if pln.propDoseCalc.useGPUtoAccumulateQuantitites
    gpu = gpuDevice();
end
switch probQuantitiesMode
    case 'phase'
        for ctIdx=1:pln.multScen.numOfCtScen
            tmpWones = ones(dij.totalNumOfBixels,1);
            tmpd = dij.physicalDoseExp{ctIdx} * tmpWones;            
            ixNz{ctIdx} = find(tmpd > 0);
        end

        matRad_cfg.dispInfo('Accumulating Omega');
        
        vCnt = 0;
        for v = voiIx
            vCnt = vCnt +1;
            matRad_cfg.dispInfo('\t\t\tStructure %d/%d',vCnt,numel(voiIx));

            for ctIdx=1:pln.multScen.numOfCtScen
                omegaCurr =  dij.physicalDoseOmega{v, ctIdx};
                structVoxels = intersect(ixNz{ctIdx},cst{v,4}{ctIdx});
                for rangeShiftScen = 1:pln.multScen.totNumRangeScen
                    if pln.multScen.scenMask(ctIdx,shiftScen,rangeShiftScen)
                        matRad_cfg.dispInfo('.');
                        scenIdx = find(ismember(pln.multScen.linearMask,[ctScen,shiftScen,rangeShiftScen], 'rows'));
                        if ~isempty(scenIdx)
                            currDij = dij.physicalDose{ctIdx,shiftScen,rangeShiftScen}(structVoxels,:);
                            omegaCurr = omegaCurr + (currDij' * currDij)*pln.multScen.scenWeight(scenIdx);
                    
                        else
                            matRad_cfg.dispWarning('Cannot accumulate quantities');
                        end
                    end
                end

                dij.physicalDoseOmega{v,ctIdx} = omegaCurr;
            end
        matRad_cfg.dispInfo('done\n');        
        end

    case 'all'
        tmpWones = ones(dij.totalNumOfBixels,1);
        tmpd = dij.physicalDoseExp{1} * tmpWones;            
        ixNz = find(tmpd > 0);
        
        matRad_cfg.dispInfo('Accumulating Omega');
        
        vCnt = 0;
        for v = voiIx
            vCnt = vCnt +1;
            matRad_cfg.dispInfo('\t\t\tStructure %d/%d',vCnt,numel(voiIx));
            structVoxels = [];
            for ctIx = pln.multScen.numOfCtScen
                newIx = intersect(ixNz,cst{v,4}{ctIx});
                structVoxels = [structVoxels; newIx];
            end
            structVoxels = unique(structVoxels);
            
            omegaCurr =  dij.physicalDoseOmega{v};
            
            for ctIdx=1:pln.multScen.numOfCtScen
                
                
                for rangeShiftScen = 1:pln.multScen.totNumRangeScen
                    if pln.multScen.scenMask(ctIdx,shiftScen,rangeShiftScen)
                        matRad_cfg.dispInfo('.');
                        scenIdx = find(ismember(pln.multScen.linearMask,[ctScen,shiftScen,rangeShiftScen], 'rows'));
                        if ~isempty(scenIdx)
                            if pln.propDoseCalc.useGPUtoAccumulateQuantitites
                                currDij = gpuArray(dij.physicalDose{ctIdx,shiftScen,rangeShiftScen}(structVoxels,:));
                                omegaCurr = omegaCurr + gather(currDij' * currDij)*(pln.multScen.scenWeight(scenIdx)./sum(pln.multScen.scenWeight));
                            else
                                currDij = dij.physicalDose{ctIdx,shiftScen,rangeShiftScen}(structVoxels,:);
                                omegaCurr = omegaCurr + (currDij' * currDij)*(pln.multScen.scenWeight(scenIdx)./sum(pln.multScen.scenWeight));
                            end

                            
                    
                        else
                            matRad_cfg.dispWarning('Cannot accumulate quantities');
                        end
                    end
                end

            end
            dij.physicalDoseOmega{v} = omegaCurr;
            matRad_cfg.dispInfo('done\n');
        end
end

dij.physicalDose(:,shiftScen,:) = cell(pln.multScen.numOfCtScen,1,pln.multScen.totNumRangeScen);

