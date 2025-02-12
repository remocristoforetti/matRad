function dij = matRad_combineMixModProbabilisticQuantities(pln,varargin)

    matRad_cfg = MatRad_Config.instance();
    
    p = inputParser();

    addParameter(p, 'protons_saveDir', [], @(x) exist(x, 'dir'));
    addParameter(p, 'photons_saveDir', [], @(x) exist(x, 'dir'));
    addParameter(p, 'carbon_saveDir',  [], @(x) exist(x, 'dir'));
    addOptional(p, 'probQuantity', []);
    
    parse(p, varargin{:});

    nModalities = pln.numOfModalities;
    modalities = {pln.originalPlans.radiationMode};
    
    if isfield(p.Results, 'probQuantity') && ~isempty(p.Results.probQuantity)
        probQuantitiesFileName = p.Results.probQuantity;
    else
        probQuantitiesFileName = 'probQuantities';
    end

    for modalityIdx=1:nModalities
    
        modalityName = modalities{modalityIdx};

        try
            currDij = load(fullfile(p.Results.([modalityName, '_saveDir']), [probQuantitiesFileName ,'.mat']), 'dij', 'expDist', 'omega');
            currDij.probQuantities.physicalDoseExp   = currDij.expDist;
            currDij.probQuantities.physicalDoseOmega = currDij.omega;
        catch
            currDij = load(fullfile(p.Results.([modalityName, '_saveDir']),  [probQuantitiesFileName ,'.mat']), 'probQuantities','dij');
        end

        fName = fieldnames(currDij.dij);
        cellFields = fName(structfun(@iscell, currDij.dij));

        emptyFields = cellfun(@(x) isempty(currDij.dij.(x){1}), cellFields);
        currDij.dij = rmfield(currDij.dij, cellFields(emptyFields));

        dij.(modalityName) = currDij.dij;
        dij.(modalityName).physicalDoseExp      = currDij.probQuantities.physicalDoseExp;

        if isfield(currDij.probQuantities, 'physicalDoseOmegaExp')        
            dij.(modalityName).physicalDoseOmegaExp = currDij.probQuantities.physicalDoseOmegaExp;
        end

        if isfield(currDij.probQuantities, 'physicalDoseJExp')
            % Check if this is a column, make it a row
            dij.(modalityName).physicalDoseJExp = cellfun(@(x) reshape(x,1, []), currDij.probQuantities.physicalDoseJExp, 'UniformOutput',false);
        end

        if isfield(currDij.probQuantities, 'alphaDoseJExp')
            dij.(modalityName).alphaDoseJExp = cellfun(@(x) reshape(x,1, []), currDij.probQuantities.alphaDoseJExp, 'UniformOutput',false);
        end

        if isfield(currDij.probQuantities, 'sqrtBetaDoseJExp')
            dij.(modalityName).sqrtBetaDoseJExp = cellfun(@(x) reshape(x,1, []), currDij.probQuantities.sqrtBetaDoseJExp, 'UniformOutput',false);
        end

        if isfield(currDij.probQuantities, 'mAlphaDoseExp')
            dij.(modalityName).mAlphaDoseExp = currDij.probQuantities.mAlphaDoseExp;
        end

        if isfield(currDij.probQuantities, 'mAlphaDoseOmegaExp')
            dij.(modalityName).mAlphaDoseOmegaExp = currDij.probQuantities.mAlphaDoseOmegaExp;
        end

        if isfield(currDij.probQuantities, 'mSqrtBetaDoseExp')
            dij.(modalityName).mSqrtBetaDoseExp = currDij.probQuantities.mSqrtBetaDoseExp;
        end

        if isfield(currDij.probQuantities, 'mSqrtBetaDoseOmegaExp')
            dij.(modalityName).mSqrtBetaDoseOmegaExp = currDij.probQuantities.mSqrtBetaDoseOmegaExp;
        end

        if isfield(currDij.probQuantities, 'alphaDoseJ')
            dij.(modalityName).alphaDoseJ = cellfun(@(x) reshape(x,1, []), currDij.probQuantities.alphaDoseJ, 'UniformOutput',false);
        end

        if isfield(currDij.probQuantities, 'mAlphaDoseOmegaCross')
            dij.(modalityName).mAlphaDoseOmegaCross = currDij.probQuantities.mAlphaDoseOmegaCross;
        end

        if isfield(currDij.probQuantities, 'mSqrtBetaDoseOmegaCross')
            dij.(modalityName).mSqrtBetaDoseOmegaCross = currDij.probQuantities.mSqrtBetaDoseOmegaCross;
        end

        if isfield(currDij.probQuantities, 'sqrtBetaDoseJ')
            dij.(modalityName).sqrtBetaDoseJ = cellfun(@(x) reshape(x,1, []), currDij.probQuantities.sqrtBetaDoseJ, 'UniformOutput',false);
        end

        if ~isfield(dij.(modalityName), 'physicalDose')
            dij.(modalityName).physicalDose = {[]};
        end
    end

    allFieldsName = cellfun(@(x) fieldnames(dij.(x)), modalities, 'UniformOutput',false);
    
    commonFieldsName = allFieldsName{1};
    for modalityIdx=2:nModalities
        commonFieldsName = intersect(commonFieldsName, allFieldsName{modalityIdx});%unique([vertcat(allFieldsName{:})]);
    end

    % Exclude large fields
    [~,excludeFieldIdx] = intersect(commonFieldsName, {'physicalDose', 'mAlphaDose', 'mSqrtBetaDose', 'mLETDose', 'physicalDoseExp', 'alphaDoseJ', 'sqrtBetaDoseJ', 'mAlphaDoseOmegaExp', 'mAlphaDoseOmegaCross'});
    commonFieldsName(excludeFieldIdx) = [];
    for propertyName = commonFieldsName'
        if isequal(dij.(modalities{1}).(propertyName{1}), dij.(modalities{2}).(propertyName{1}))
            dij.(propertyName{1}) = dij.(modalities{1}).(propertyName{1});
        else
            dij.(propertyName{1}) = cellfun(@(x) dij.(x).(propertyName{1}),modalities, 'UniformOutput', false);
        end
    end
     
    dij.numOfModalities = nModalities;
    dij.radiationModalities = modalities;
    
    dij.totalNumOfBixels = sum([dij.totalNumOfBixels{:}]);
end