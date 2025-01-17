function [physicalDose, physicalDoseOmegaReduced, mAlphaDose, mSqrtBetaDose, mLETDose,mLETdJ, alphaJ, sqrtBetaJ,physicalDoseJ]  = matRad_loadScenariosFromMeta(saveDir, scenariosMeta, dijTemplate, verbosity)

    matRad_cfg = MatRad_Config.instance();
    %% Input
    if ~exist('dijTemplate', 'var') || isempty(dijTemplate)

        matRad_cfg.dispInfo('No dijTemplate provided, trying to locate one ...');
        try
            load(fullfile(saveDirectory, 'dijTemplate.mat'),'dijTemplate');
            matRad_cfg.dispInfo('done. \n');
        catch

            matRad_cfg.dispError('Unable to load dijTemplate file');
        end
    end

    if ~exist('verbosity', 'var') || isempty(verbosity)

        verbosity = false;
    end

    %% Output
    %nOutputs = nargout();
    %TODO
    %% Loading
    physicalDose             = {};
    physicalDoseOmegaReduced = {};
    mAlphaDose               = {};
    mSqrtBetaDose            = {};
    mLETDose                 = {};
    mLETdJ                   = {};
    alphaJ                   = {};
    sqrtBetaJ                = {};
    physicalDoseJ            = {};

    % Get meta info from dij template
    nVoxels = dijTemplate.doseGrid.numOfVoxels;
    nBixels = dijTemplate.totalNumOfBixels;

    nScensToLoad = numel(scenariosMeta);
    
    physicalDose = arrayfun(@(scen) spalloc(nVoxels,nBixels,scen.nnzElements), scenariosMeta, 'UniformOutput',false);
    %mAlphaDose = arrayfun(@(scen) spalloc(nVoxels,nBixels,scen.nnzElements), scenariosMeta, 'UniformOutput',false);
    %mSqrtBetaDose = arrayfun(@(scen) spalloc(nVoxels,nBixels,scen.nnzElements), scenariosMeta, 'UniformOutput',false);

    alphaJ        = cell(nScensToLoad,1);
    sqrtBetaJ     = cell(nScensToLoad,1);
    % sqrtBetaJ     = arrayfun(@(scen) spalloc(nBixels,1, nBixels), scenariosMeta, 'UniformOutput',false);
    
    
    stringLength = 0;
    for scenIdx=1:nScensToLoad

        if verbosity
            fprintf(repmat('\b',1,stringLength));
            stringLength = fprintf('\tLoading scenario: %u/%u\n', scenIdx, nScensToLoad);
        end
        
        currScenarioMeta = scenariosMeta(scenIdx);
        fileName = fullfile(saveDir, currScenarioMeta.name);
        
        currDijScen = load(fileName, 'dijScenario');
        currDijScen = currDijScen.dijScenario;

        if ~isstruct(currDijScen)
            % This is for older compatibility
            physicalDose(scenIdx) = currDijScen;
        else
            physicalDose(scenIdx)  = currDijScen.physicalDose;
            
            if isfield(currDijScen, 'mAlphaDose')
                mAlphaDose{scenIdx}   = currDijScen.mAlphaDose;
            end

            if isfield(currDijScen, 'mSqrtBetaDose')
                mSqrtBetaDose{scenIdx} = currDijScen.mSqrtBetaDose;
            end

            if isfield(currDijScen, 'mLETDose')
                mLETDose{scenIdx} = currDijScen.mLETDose;
            end

            if isfield(currDijScen, 'alphaJ') && isfield(currDijScen, 'sqrtBetaJ')
                alphaJ{scenIdx}    = currDijScen.alphaJ;
                sqrtBetaJ{scenIdx} = currDijScen.sqrtBetaJ;
            end
            
            if isfield(currDijScen, 'physicalDoseOmegaReduced')
                physicalDoseOmegaReduced{scenIdx} = currDijScen.physicalDoseOmegaReduced;
            else
                physicalDoseOmegaReduced{scenIdx} = [];
            end

            if isfield(currDijScen, 'mLETdJ')
                mLETdJ{scenIdx} = currDijScen.mLETdJ;
            end

            
            if isfield(currDijScen, 'mLETdJ')
                physicalDoseJ{scenIdx} = currDijScen.physicalDoseJ;
            end

        end
    end
end