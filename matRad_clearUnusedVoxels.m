function [dij, originalDij] = matRad_clearUnusedVoxels(cst, dij)

matRad_cfg = MatRad_Config.instance();

matRad_cfg.dispInfo('Disabling Voxels outside of ROIs');


%check for structures that do have objective/constraint
includeMask = zeros(dij.doseGrid.dimensions);

for  i = 1:size(cst,1)
    
    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )
        
        % loop over the number of constraints for the current VOI
        for j = 1:numel(cst{i,6})
            
            objective = cst{i,6}{j};
            
            % only perform gradient computations for objectiveectives
            if isa(objective,'DoseObjectives.matRad_DoseObjective') || isa(objective, 'DoseConstraints.matRad_DoseConstraint')

                includeMask(cst{i,4}{1}) = 1;

            end

        end

    end

end

useScen = [1:numel(dij.physicalDose)];

% v = zeros(size(dij.physicalDose{1},1),1);
% w = ones(size(dij.physicalDose{1},2),1);
% tic
% for k=1:10
%     v = dij.physicalDose{1}*w;
% end
% originalTime = toc;


if isfield(dij, 'physicalDose')
    for scenIdx = useScen
        % Preserve original dij for calcCubes;
        originalDij.physicalDose{scenIdx} = dij.physicalDose{scenIdx};
        dij.physicalDose{scenIdx} = dij.physicalDose{scenIdx}.*includeMask(:);
    end
end

if isfield(dij, 'mAlphaDose') && isfield(dij, 'mSqrtBetaDose')
    for scenIdx = useScen

        originalDij.mAlphaDose{scenIdx} = dij.mAlphaDose{scenIdx};
        originalDij.mSqrtBetaDose{scenIdx} = dij.mSqrtBetaDose{scenIdx};

        dij.mAlphaDose{scenIdx} = dij.mAlphaDose{scenIdx}.*includeMask(:);
        dij.mSqrtBetaDose{scenIdx} = dij.mSqrtBetaDose{scenIdx}.*includeMask(:);
    end
end

originalSparsity = arrayfun(@(scen) nnz(originalDij.physicalDose{scen}), useScen, 'UniformOutput',false);
reducedDijSparsity = arrayfun(@(scen) nnz(dij.physicalDose{scen}), useScen, 'UniformOutput',false);

matRad_cfg.dispInfo('... done\n');

if  ~any([originalSparsity{:}] > [reducedDijSparsity{:}])
    matRad_cfg.dispWarning('No dij sparsity reduction obtained, there is no benefit in reducing the dij volume');
end
% w = ones(size(dij.physicalDose{1},2),1);
% tic
% for k=1:10
%     v = dij.physicalDose{1}*w;
% end
% afterTime = toc;

end