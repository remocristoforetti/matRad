function dij = matRad_combineMMdijs(pln, varargin)
    
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
    
        currDij = load(fullfile(p.Results.([modalityName, '_saveDir']),  [probQuantitiesFileName ,'.mat']), 'dij');
    
        currNominalScenMeta = load(fullfile(p.Results.([modalityName, '_saveDir']), 'scenario_1.mat'), 'absRangeShift', 'isoShift', 'relRangeShift');
    
        if all([currNominalScenMeta.absRangeShift == 0, currNominalScenMeta.isoShift == [0,0,0], currNominalScenMeta.relRangeShift == 0])
            currNominalScenario = load(fullfile(p.Results.([modalityName, '_saveDir']), 'scenario_1.mat'), 'dijScenario');
        else
            matRad_cfg.error('Unable to find nominal scenario');
        end
    
        fName = fieldnames(currDij.dij);
        cellFields = fName(structfun(@iscell, currDij.dij));
    
        emptyFields = cellfun(@(x) isempty(currDij.dij.(x){1}), cellFields);
        currDij.dij = rmfield(currDij.dij, cellFields(emptyFields));
    
        dij.(modalityName) = currDij.dij;
    
        scenFieldNames = fieldnames(currNominalScenario.dijScenario);
        for scenFName=scenFieldNames'
            %if isfield(currNominalScenario.dijScenario, scenFName{1})
                dij.(modalityName).(scenFName{1}) = currNominalScenario.dijScenario.(scenFName{1});
            %end
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