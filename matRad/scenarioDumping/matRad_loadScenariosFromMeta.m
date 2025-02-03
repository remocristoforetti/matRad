function [outScenarioQuantities]  = matRad_loadScenariosFromMeta(saveDir, scenariosMeta, dijTemplate, bioQuantities, verbosity)

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

    if ~exist('bioQuantities', 'var') || isempty(bioQuantities)

        bioQuantities = false;
    end

    if ~exist('verbosity', 'var') || isempty(verbosity)
        verbosity = false;
    end

    %% Loading
    physicalDose             = {};

    if bioQuantities
        mAlphaDose               = {};
        mSqrtBetaDose            = {};
        alphaJ                   = {};
        sqrtBetaJ                = {};
        physicalDoseJ            = {};
    end

    mLETDose                 = {};
    mLETdJ                   = {};

    % Get meta info from dij template
    nVoxels = dijTemplate.doseGrid.numOfVoxels;
    nBixels = dijTemplate.totalNumOfBixels;

    nScensToLoad = numel(scenariosMeta);
    
    physicalDose = arrayfun(@(scen) spalloc(nVoxels,nBixels,scen.nnzElements), scenariosMeta, 'UniformOutput',false);
    if bioQuantities
        mAlphaDose = arrayfun(@(scen) spalloc(nVoxels,nBixels,scen.nnzElements), scenariosMeta, 'UniformOutput',false);
        mSqrtBetaDose = arrayfun(@(scen) spalloc(nVoxels,nBixels,scen.nnzElements), scenariosMeta, 'UniformOutput',false);
    
        alphaJ        = cell(nScensToLoad,1);
        sqrtBetaJ     = cell(nScensToLoad,1);
    end
    
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
            
            if bioQuantities
                if isfield(currDijScen, 'mAlphaDose')
                    mAlphaDose(scenIdx)   = currDijScen.mAlphaDose;
                end

                if isfield(currDijScen, 'mSqrtBetaDose')
                    mSqrtBetaDose(scenIdx) = currDijScen.mSqrtBetaDose;
                end
                
                % if isfield(currDijScen, 'alphaDoseJ') && isfield(currDijScen, 'sqrtBetaDoseJ')
                %     alphaJ(scenIdx)    = currDijScen.alphaDoseJ;
                %     sqrtBetaJ(scenIdx) = currDijScen.sqrtBetaDoseJ;
                % end

            end

            if isfield(currDijScen, 'mLETDose')
                mLETDose(scenIdx) = currDijScen.mLETDose;
            end


            % if isfield(currDijScen, 'mLETdJ')
            %     mLETdJ(scenIdx) = currDijScen.mLETdJ;
            % end

            
            % if isfield(currDijScen, 'mLETdJ')
            %     physicalDoseJ(scenIdx) = currDijScen.physicalDoseJ;
            % end

        end
    end

    outScenarioQuantities.physicalDose  = physicalDose;
    % outScenarioQuantities.physicalDoseJ = physicalDoseJ;

    if bioQuantities
        outScenarioQuantities.mAlphaDose = mAlphaDose;
        outScenarioQuantities.mSqrtBetaDose = mSqrtBetaDose;
        % outScenarioQuantities.alphaDoseJ = alphaJ;
        % outScenarioQuantities.sqrtBetaDoseJ = sqrtBetaJ;
    end


    outScenarioQuantities.mLETdose = mLETDose;
    % outScenarioQuantities.mLETdJ   = mLETdJ;
end