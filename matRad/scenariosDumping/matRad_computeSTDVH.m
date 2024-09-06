function [stdvh, expDist, stdDist, stdGrid] = matRad_computeSTDVH(cst, scenariosDistributions, weights, stdGrid, fractions)
%   matRad_computeSTDVH(cst, scenariosDistributions, weights, stdGrid, fractions)
%   scenarioDistributions = cella raay with dose distributions
%


    matRad_cfg = MatRad_Config.instance();
    if ~exist('weights', 'var')
        weights = ones(numel(scenariosDistributions),1);
    end
    
    if ~iscolumn(weights)
        weights = weights';
    end
    
    if ~exist('fractions', 'var')
        fractions = 1;
    end

    originalDimensions = size(scenariosDistributions{1});

    for scenIdx=1:numel(scenariosDistributions)
        allScen(:,scenIdx) = scenariosDistributions{scenIdx}(:)*fractions;
    end
    
    expDist = allScen*weights;
    
    stdDist = std(allScen, weights,2);

    expDist = reshape(expDist, originalDimensions);
    stdDist = reshape(stdDist, originalDimensions);
    
    if ~exist('stdGrid', 'var') || isempty(stdGrid)
        stdGrid = linspace(0,1.1*max(stdDist, [], 'all'),1000);

    end

    matRad_cfg.dispWarning('SDVH estimated for reference CT scenario #%i', 1);
    stdvh = matRad_calcDVH(cst,stdDist,1,[],stdGrid);
end