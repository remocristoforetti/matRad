function [includeMask] = matRad_getRobustVoxelsOnGrid(cst, doseGrid, VdoseGrid)
    matRad_cfg = MatRad_Config.instance();
    
    matRad_cfg.dispInfo('Disabling dose calculation in voxels outside of ROIs in robustness scenarios');
    cst  = matRad_setOverlapPriorities(cst);
    
    %check for structures that do have objective/constraint
    includeMask = cell(size(cst{1,4},2),1);
    includeMask(:) = {zeros(prod(doseGrid.dimensions),1)};
    %maskT = zeros(prod(doseGrid.dimensions),1);
    for  i = 1:size(cst,1)
        
            % Only take OAR or target VOI.
            if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )
                
                % loop over the number of constraints for the current VOI
                for j = 1:numel(cst{i,6})
                    
                    obj = cst{i,6}{j};
                     
                    if ~isa(obj,'matRad_DoseOptimizationFunction')
                        try
                            obj = matRad_DoseOptimizationFunction.createInstanceFromStruct(obj);
                        catch
                            matRad_cfg.dispError('cst{%d,6}{%d} is not a valid Objective/constraint! Remove or Replace and try again!',i,j);
                        end
                    end
                    
                    
                    robustness = obj.robustness;
                    %maskT(cst{i,4}{ctIdx}) = i;
                    if isa(obj,'DoseObjectives.matRad_DoseObjective') || isa(obj, 'DoseConstraints.matRad_DoseConstraint')
                        

                        if ~isequal(robustness, 'none')
                            for ctIdx=1:size(cst{i,4},2)
                                includeMask{ctIdx}(cst{i,4}{ctIdx}) = 1;
                                
                            end
                        else
                            matRad_cfg.dispInfo('\n \t \t no robustness specification detected in cst structure %s.', cst{i,2});
                        end
                    end
                    
        
                end
        
            end
 
    end

    for ctIdx=1:size(cst{i,4},2)
        includeMask{ctIdx} = includeMask{ctIdx}(VdoseGrid);
    end
    matRad_cfg.dispInfo('\n done\n');
end