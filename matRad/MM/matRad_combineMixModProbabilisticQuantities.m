function dij = matRad_combineMixModProbabilisticQuantities(pln,varargin)

    matRad_cfg = MatRad_Config.instance();
    
    p = inputParser();

    addParameter(p, 'protons_saveDir', [], @(x) exist(x, 'dir'));
    addParameter(p, 'photons_saveDir', [], @(x) exist(x, 'dir'));
    addParameter(p, 'carbon_saveDir',  [], @(x) exist(x, 'dir'));

    parse(p, varargin{:});

    nModalities = pln.numOfModalities;
    modalities = {pln.originalPlans.radiationMode};
    

    for modalityIdx=1:nModalities
    
        modalityName = modalities{modalityIdx};

        try
            currDij = load(fullfile(p.Results.([modalityName, '_saveDir']), 'probQuantities.mat'), 'dij', 'expDist', 'omega');
            currDij.probQuantities.physicalDoseExp   = currDij.expDist;
            currDij.probQuantities.physicalDoseOmega = currDij.omega;
        catch
            currDij = load(fullfile(p.Results.([modalityName, '_saveDir']), 'probQuantities.mat'), 'probQuantities','dij');
        end

        fName = fieldnames(currDij.dij);
        cellFields = fName(structfun(@iscell, currDij.dij));

        emptyFields = cellfun(@(x) isempty(currDij.dij.(x){1}), cellFields);
        currDij.dij = rmfield(currDij.dij, cellFields(emptyFields));

        dij.(modalityName) = currDij.dij;
        dij.(modalityName).physicalDose   = currDij.probQuantities.physicalDoseExp;
        dij.(modalityName).physicalDoseOmega = currDij.probQuantities.physicalDoseOmega;

        if isfield(currDij.probQuantities, 'mAlphaDoseExp')
            dij.(modalityName).mAlphaDose = currDij.probQuantities.mAlphaDoseExp;
        end

        if isfield(currDij.probQuantities, 'mAlphaDoseOmega')
            dij.(modalityName).mAlphaDoseOmega = currDij.probQuantities.mAlphaDoseOmega;
        end

        if isfield(currDij.probQuantities, 'mSqrtBetaDoseExp')
            dij.(modalityName).mSqrtBetaDose = currDij.probQuantities.mSqrtBetaDoseExp;
        end

        if isfield(currDij.probQuantities, 'mSqrtBetaDoseOmega')
            dij.(modalityName).mSqrtBetaDoseOmega = currDij.probQuantities.mSqrtBetaDoseOmega;
        end

        if isfield(currDij.probQuantities, 'alphaDoseJ')
            dij.(modalityName).alphaDoseJ = currDij.probQuantities.alphaDoseJ;
        end

        if isfield(currDij.probQuantities, 'sqrtBetaDoseJ')
            dij.(modalityName).sqrtBetaDoseJ = currDij.probQuantities.sqrtBetaDoseJ;
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
    [~,excludeFieldIdx] = intersect(commonFieldsName, {'physicalDose', 'mAlphaDose', 'mSqrtBetaDose', 'mLETDose', 'physicalDoseExp', 'alphaDoseJ', 'sqrtBetaDoseJ'});
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