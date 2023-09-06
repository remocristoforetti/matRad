function [includeMask] = matRad_getRobustVoxelsOnGrid(cst, doseGrid, VdoseGrid, selectionMode)
    %voxels set to 1 in includeMask will have dose calculation

    matRad_cfg = MatRad_Config.instance();
    
    matRad_cfg.dispInfo('Disabling dose calculation in voxels outside of ROIs in robustness scenarios');
    cst  = matRad_setOverlapPriorities(cst);

    selectedCstStructs = [];
    
    includeMask = cell(size(cst{1,4},2),1);
    includeMask(:) = {zeros(prod(doseGrid.dimensions),1)};
        
    %structures that are selected here will be included in dose calculation
    %over the robust scenarios
    if isequal(selectionMode , 'none')
        for ctScenIdx=1:size(includeMask,2)
            includeMask{ctScenIdx}(:) = 1;
        end
    else

        switch selectionMode
        
            case 'all'
                selectedCstStructs = [1:size(cst,1)];
            case 'targetOnly'
                selectedCstStructs = find(cellfun(@(x) strcmp(x,'TARGET'), [cst(:,3)]));
            case 'objectivesOnly'
                for i=1:size(cst,1)
                    if numel(cst{i,6})>0
                        selectedCstStructs = [selectedCstStructs, i];
                    end
                end
            case 'oarsOnly'
                selectedCstStructs = find(cellfun(@(x) strcmp(x,'OAR'), [cst(:,3)]));    
            otherwise
        
        end
    
        
        %check for structures that do have objective/constraint
        
    
    
        %loop over all cst sturctures 
        for i=1:size(cst,1)
        
            if ~isempty(cst{i,4}{1})
                
                if numel(cst{i,6}) > 0
                    %loop over obj/constraint functions
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
    
                        if any(intersect(i, selectedCstStructs))
                            for ctIdx=1:size(cst{i,4},2)
                                includeMask{ctIdx}(cst{i,4}{ctIdx}) = 1;
                            end
    
                            if isequal(robustness, 'none')
                                matRad_cfg.dispWarning('Dose calculation is performed on cst structure %s, but this structure has no robustness.', cst{i,2});
                            end
                        else
                            if ~isequal(robustness, 'none')
                                matRad_cfg.dispError('Trying to clear voxels in cst structure: %s, but this structure has robustness requirements', cst{i,2});
                            end
                        end
                    end
    
                else %numel(cst{i,6}) < 0
                    if any(intersect(i, selectedCstStructs))
                        matRad_cfg.dispWarning('Trying to clear voxels in cst structure: %s, but this structure does not have any objective or constraint', cst{i,2}');
                    end
                end %numel(cst{i,6}) > 0
            end %if cst{i,4} not empty
        
        end %for loop over cst
    
    
        % for  i = selectedCstStructs
        % 
        %         % Only take OAR or target VOI.
        %         if ~isempty(cst{i,4}{1})
        % 
        %             % loop over the number of constraints for the current VOI
        %             if numel(cst{i,6}) > 0
        %                 for j = 1:numel(cst{i,6})
        % 
        %                     obj = cst{i,6}{j};
        % 
        %                     if ~isa(obj,'matRad_DoseOptimizationFunction')
        %                         try
        %                             obj = matRad_DoseOptimizationFunction.createInstanceFromStruct(obj);
        %                         catch
        %                             matRad_cfg.dispError('cst{%d,6}{%d} is not a valid Objective/constraint! Remove or Replace and try again!',i,j);
        %                         end
        %                     end
        % 
        % 
        %                     robustness = obj.robustness;
        % 
        %                     if ~isequal(robustness, 'none')
        %                         for ctIdx=1:size(cst{i,4},2)
        %                             includeMask{ctIdx}(cst{i,4}{ctIdx}) = 1;
        % 
        %                         end
        %                     else
        %                         matRad_cfg.dispInfo('\n \t \t no robustness specification detected in cst structure %s.', cst{i,2});
        %                     end
        % 
        % 
        %                 end
        %             else
        %                 matRad_cfg.dispWarning('Trying to clean voxels from structure %s, but no optimization function is associated to this structure', cst{i,2});
        %             end
        % 
        %         end
        % 
        % end
    end

    for ctIdx=1:size(cst{1,4},2)
        includeMask{ctIdx} = includeMask{ctIdx}(VdoseGrid);
    end
    
    matRad_cfg.dispInfo('\n done\n');
end