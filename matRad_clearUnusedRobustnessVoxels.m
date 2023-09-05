function [dij, originalDij] = matRad_clearUnusedRobustnessVoxels(cst, dij, robustMask)


% Exclude 
useScen = [2:numel(dij.physicalDose)];

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