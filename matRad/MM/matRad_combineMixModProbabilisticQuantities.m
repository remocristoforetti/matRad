function dij = matRad_combineMixModProbabilisticQuantities(pln,varargin)

    matRad_cfg = MatRad_Config.instance();
    
    defaultFieldsToLoad = {'physicalDoseOmegaExp',... % all
                            'physicalDoseJExp',...
                            'alphaDoseJExp',...
                            'sqrtBetaDoseJExp',...
                            'mAlphaDoseExp',...
                            'mAlphaDoseOmegaExp',...
                            'mSqrtBetaDoseExp',...
                            'mSqrtBetaDoseOmegaExp',...
                            'alphaDoseJ',...
                            'mAlphaDoseOmegaCross',...
                            'mSqrtBetaDoseOmegaCross',...
                            'sqrtBetaDoseJ',...
                            'mEffectOmegaHP'};
    p = inputParser();

    addParameter(p, 'protons_saveDir', [], @(x) exist(x, 'dir'));
    addParameter(p, 'photons_saveDir', [], @(x) exist(x, 'dir'));
    addParameter(p, 'carbon_saveDir',  [], @(x) exist(x, 'dir'));
    addOptional(p, 'probQuantity', []);
    addOptional(p, 'fieldsToLoad', defaultFieldsToLoad);
    
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

        % Load filedss into the dij
        fieldsToLoad = p.Results.fieldsToLoad;
        for i = 1:length(fieldsToLoad)
            field = fieldsToLoad{i};
            if isfield(currDij.probQuantities, field)
                if iscell(currDij.probQuantities.(field))
                    % dij.(modalityName).(field) = cellfun(@(x) reshape(x, 1, []), currDij.probQuantities.(field), 'UniformOutput', false);
                    dij.(modalityName).(field) = currDij.probQuantities.(field);%cellfun(@(x) reshape(x, 1, []), currDij.probQuantities.(field), 'UniformOutput', false);
                    currDij.probQuantities.(field) = [];
                else
                    dij.(modalityName).(field) = currDij.probQuantities.(field);
                end
            end
        end
    end

    allFieldsName = cellfun(@(x) fieldnames(dij.(x)), modalities, 'UniformOutput',false);
    
    commonFieldsName = allFieldsName{1};
    for modalityIdx=2:nModalities
        commonFieldsName = intersect(commonFieldsName, allFieldsName{modalityIdx});
    end

    % Exclude large fields
    [~,excludeFieldIdx] = intersect(commonFieldsName, {'physicalDose', 'mAlphaDose', 'mSqrtBetaDose', 'mLETDose', 'physicalDoseExp', 'alphaDoseJ', 'sqrtBetaDoseJ', 'mAlphaDoseOmegaExp', 'mAlphaDoseOmegaCross', 'mEffectOmegaHP'});
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