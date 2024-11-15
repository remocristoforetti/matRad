function plnJO = matRad_plnWrapper(pln)
% Combines the arbitrary number of input plans into a single plan for the
% different modalities.
% 
% call
%   plnJO = matRad_plnWrapper(pln)
%
% input
%   pln:       array of pln structure for the different modalities (if any)
%
% output
%   plnJO:      synthetic overarching pln stuct for Joint Opt 
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matRad_cfg = MatRad_Config.instance();
nPlans = length(pln);
if nPlans>1
    
    % Collect all fields from the single modelity plans
    allFields = arrayfun(@(plan) fieldnames(plan), pln, 'UniformOutput',false);
    allFields = unique(cat(1,allFields{:}));
    plnJO = cell2struct(cell(1,length(allFields)), allFields, 2);
   
    % No need for plan consistency anymore with new biological models, will
    % just check later
    %pln = matRad_plnConsistency(pln); %to be reviewed
   
    % Save now the plns. This will have now the updated pln.bioParam
    originalPlans = pln;
   
    %Define the Pln properties
    plnJO.numOfFractions  = sum([pln(:).numOfFractions]); %Use total number of fractions
    plnJO.radiationMode   = 'MixMod';
    plnJO.machine         = 'MixMod';
    plnJO.propStf         = [pln(:).propStf];
    plnJO.numOfModalities = nPlans;

    %%%%%%% propDoseCalc %%%%%
    allFields = arrayfun(@(plan) fieldnames(plan.propDoseCalc), pln, 'UniformOutput',false);
    allFields = unique(cat(1,allFields{:}));
    plnJO.propDoseCalc(1:plnJO.numOfModalities) = cell2struct(cell(1,length(allFields)), allFields, 2);

    for modalityIdx=1:plnJO.numOfModalities
        currModalityFields = fieldnames(pln(modalityIdx).propDoseCalc);
        for fieldIdx=1:length(currModalityFields)
           plnJO.propDoseCalc(modalityIdx).(currModalityFields{fieldIdx}) = pln(modalityIdx).propDoseCalc.(currModalityFields{fieldIdx});
        end
    end

    %%%%%%% propStf %%%%%
    allFields = arrayfun(@(plan) fieldnames(plan.propStf), pln, 'UniformOutput',false);
    allFields = unique(cat(1,allFields{:}));
    plnJO.propStf(1:plnJO.numOfModalities) = cell2struct(cell(1,length(allFields)), allFields, 2);

    for modalityIdx=1:plnJO.numOfModalities
        currModalityFields = fieldnames(pln(modalityIdx).propStf);
        for fieldIdx=1:length(currModalityFields)
           plnJO.propStf(modalityIdx).(currModalityFields{fieldIdx}) = pln(modalityIdx).propStf.(currModalityFields{fieldIdx});
        end
    end

    %%%%%%% propOpt %%%%%
    allFields = arrayfun(@(plan) fieldnames(plan.propOpt), pln, 'UniformOutput',false);
    allFields = unique(cat(1,allFields{:}));
    plnJO.propOpt = cell2struct(cell(1,length(allFields)), allFields, 2);
    
    
    for modalityIdx=1:plnJO.numOfModalities

        currModalityFields = fieldnames(pln(modalityIdx).propOpt);
        for fieldIdx=1:length(currModalityFields)
            if ~any(strcmp(currModalityFields{fieldIdx}, {'spatioTemp', 'STfractions', 'STscenarios'}))
                plnJO.propOpt.(currModalityFields{fieldIdx}) = pln(modalityIdx).propOpt.(currModalityFields{fieldIdx});
            end
        end
    end

    %%% Disable STfractionation fo the time being
    for modalityIdx=1:plnJO.numOfModalities
        if (isfield(pln(modalityIdx).propOpt, 'spatioTemp')) && pln(modalityIdx).propOpt.spatioTemp
            matRad.dispWarning('Sorry no spatioTemporal avalability yet');
        end

    end
    plnJO.propOpt.spatioTemp = zeros(1,plnJO.numOfModalities);
    plnJO.propOpt.STfractions = {pln.numOfFractions};
    plnJO.propOpt.STscenarios = ones(1,plnJO.numOfModalities);

    % Feed the first bio model quantity, they are consistent
    plnJO.bioParam = matRad_bioModel(plnJO.radiationMode, 'MixMod');
    plnJO.bioParam.singleModalityModels = {pln(:).bioParam};
    plnJO.originalPlans = originalPlans;
    plnJO.multScen = [pln(:).multScen];

else
   %Do nothing
   plnJO = pln;
end


end