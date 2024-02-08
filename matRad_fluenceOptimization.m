function [resultGUI,optimizer] = matRad_fluenceOptimization(dij,cst,pln,wInit)
% matRad inverse planning wrapper function
%
% call
%   [resultGUI,optimizer] = matRad_fluenceOptimization(dij,cst,pln)
%   [resultGUI,optimizer] = matRad_fluenceOptimization(dij,cst,pln,wInit)
%
% input
%   dij:        matRad dij struct
%   cst:        matRad cst struct
%   pln:        matRad pln struct
%   wInit:      (optional) custom weights to initialize problems
%
% output
%   resultGUI:  struct containing optimized fluence vector, dose, and (for
%               biological optimization) RBE-weighted dose etc.
%   optimizer:  Used Optimizer Object
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

% consider VOI priorities
cst  = matRad_setOverlapPriorities(cst);

% check & adjust objectives and constraints internally for fractionation
for i = 1:size(cst,1)
    %Compatibility Layer for old objective format
    if isstruct(cst{i,6})
        cst{i,6} = arrayfun(@matRad_DoseOptimizationFunction.convertOldOptimizationStruct,cst{i,6},'UniformOutput',false);
    end
    for j = 1:numel(cst{i,6})

        obj = cst{i,6}{j};

        %In case it is a default saved struct, convert to object
        %Also intrinsically checks that we have a valid optimization
        %objective or constraint function in the end
        if ~isa(obj,'matRad_DoseOptimizationFunction') && ~isa(obj,'OmegaObjectives.matRad_OmegaObjective')
            try
                obj = matRad_DoseOptimizationFunction.createInstanceFromStruct(obj);
            catch
                matRad_cfg.dispError('cst{%d,6}{%d} is not a valid Objective/constraint! Remove or Replace and try again!',i,j);
            end
        end

        if isa(obj, 'matRad_DoseOptimizationFunction')
            obj = obj.setDoseParameters(obj.getDoseParameters()/pln.numOfFractions);
        end
        cst{i,6}{j} = obj;
    end

end

% resizing cst to dose cube resolution
cst = matRad_resizeCstToGrid(cst,dij.ctGrid.x,  dij.ctGrid.y,  dij.ctGrid.z,...
    dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);

% find target indices and described dose(s) for weight vector
% initialization
V          = [];
doseTarget = [];
ixTarget   = [];

for i = 1:size(cst,1)
    if isequal(cst{i,3},'TARGET') && ~isempty(cst{i,6})
        V = [V;cst{i,4}{1}];

        %Iterate through objectives/constraints
        fDoses = [];
        for fObjCell = cst{i,6}
            dParams = fObjCell{1}.getDoseParameters();
            %Don't care for Inf constraints
            dParams = dParams(isfinite(dParams));
            %Add do dose list
            fDoses = [fDoses dParams];
        end

        doseTarget = [doseTarget fDoses];
        ixTarget   = [ixTarget i*ones(1,length(fDoses))];
    end

end
[doseTarget,i] = max(doseTarget);
ixTarget       = ixTarget(i);
wOnes          = ones(dij.totalNumOfBixels,1);

% calculate initial beam intensities wInit
matRad_cfg.dispInfo('Estimating initial weights... ');

if exist('wInit','var')
    %do nothing as wInit was passed to the function
    matRad_cfg.dispInfo('chosen provided wInit!\n');

    % Write ixDose which is needed for the optimizer
    if pln.bioParam.bioOpt
        dij.ixDose  = dij.bx~=0;

        %pre-calculations
        dij.gamma             = zeros(dij.doseGrid.numOfVoxels,dij.numOfScenarios);
        dij.gamma(dij.ixDose) = dij.ax(dij.ixDose)./(2*dij.bx(dij.ixDose));
    end


elseif strcmp(pln.bioParam.model,'constRBE') && strcmp(pln.radiationMode,'protons')
    % check if a constant RBE is defined - if not use 1.1
    if ~isfield(dij,'RBE')
        dij.RBE = 1.1;
    end

    if ~isempty(dij.physicalDose{1})
        doseTmp = dij.physicalDose{1}*wOnes;
    else
        doseTmp = dij.physicalDoseExp{1}*wOnes;
    end
    bixelWeight =  (doseTarget)/(dij.RBE * mean(doseTmp(V)));
    wInit       = wOnes * bixelWeight;
    matRad_cfg.dispInfo('chosen uniform weight of %f!\n',bixelWeight);

elseif pln.bioParam.bioOpt
    % retrieve photon LQM parameter
    [ax,bx] = matRad_getPhotonLQMParameters(cst,dij.doseGrid.numOfVoxels,dij.numOfScenarios);

    if ~isequal(dij.ax(dij.ax~=0),ax(dij.ax~=0)) || ...
            ~isequal(dij.bx(dij.bx~=0),bx(dij.bx~=0))
        matRad_cfg.dispError('Inconsistent biological parameter - please recalculate dose influence matrix!\n');
    end

    for i = 1:size(cst,1)

        for j = 1:size(cst{i,6},2)
            % check if prescribed doses are in a valid domain
            if any(cst{i,6}{j}.getDoseParameters() > 5) && isequal(cst{i,3},'TARGET')
                matRad_cfg.dispError('Reference dose > 5 Gy[RBE] for target. Biological optimization outside the valid domain of the base data. Reduce dose prescription or use more fractions.\n');
            end

        end
    end

    dij.ixDose  = dij.bx~=0;

    if isequal(pln.bioParam.quantityOpt,'effect')

        effectTarget = cst{ixTarget,5}.alphaX * doseTarget + cst{ixTarget,5}.betaX * doseTarget^2;
        aTmp = dij.mAlphaDose{1}*wOnes;
        bTmp = dij.mSqrtBetaDose{1} * wOnes;
        p = sum(aTmp(V)) / sum(bTmp(V).^2);
        q = -(effectTarget * length(V)) / sum(bTmp(V).^2);

        wInit        = -(p/2) + sqrt((p^2)/4 -q) * wOnes;

    elseif isequal(pln.bioParam.quantityOpt,'RBExD')

        %pre-calculations
        dij.gamma             = zeros(dij.doseGrid.numOfVoxels,dij.numOfScenarios);
        dij.gamma(dij.ixDose) = dij.ax(dij.ixDose)./(2*dij.bx(dij.ixDose));

        % calculate current effect in target
        aTmp = dij.mAlphaDose{1}*wOnes;
        bTmp = dij.mSqrtBetaDose{1} * wOnes;
        doseTmp = dij.physicalDose{1}*wOnes;

        CurrEffectTarget = aTmp(V) + bTmp(V).^2;
        % ensure a underestimated biological effective dose
        TolEstBio        = 1.2;
        % calculate maximal RBE in target
        maxCurrRBE = max(-cst{ixTarget,5}.alphaX + sqrt(cst{ixTarget,5}.alphaX^2 + ...
            4*cst{ixTarget,5}.betaX.*CurrEffectTarget)./(2*cst{ixTarget,5}.betaX*doseTmp(V)));
        wInit    =  ((doseTarget)/(TolEstBio*maxCurrRBE*max(doseTmp(V))))* wOnes;
    end

    matRad_cfg.dispInfo('chosen weights adapted to biological dose calculation!\n');
else
    if ~isempty(dij.physicalDose{1})
        doseTmp = dij.physicalDose{1}*wOnes;
    else
        doseTmp = dij.physicalDoseExp{1}*wOnes;
    end
    bixelWeight =  (doseTarget)/mean(doseTmp(V));
    wInit       = wOnes * bixelWeight;
    matRad_cfg.dispInfo('chosen uniform weight of %f!\n',bixelWeight);
end


%% calculate probabilistic quantities for probabilistic optimization if at least

%Get all non-empty scenario indexes from physicalDose
allScen = find(~cellfun(@isempty, dij.physicalDose));

%If no 4D optimization, scen4D is set to 1
if isfield(pln,'propOpt') && isfield(pln.propOpt,'scen4D')
    scen4D = pln.propOpt.scen4D;
else
    scen4D = 1;
end

if isequal(scen4D,'all')
    scen4D = [1:pln.multScen.numOfCtScen];
end

%Set here useScen to all thze nominal ct scenarios selected, if there is
%robustness or no need for scenario calculation, useScen will be changed
%accordingly
nominalCTScen = scen4D;
useScen = nominalCTScen;

%Check for robust optimization and probabilistic optimization requirements
PROB_FLAG = false;
ROB_FLAG = false;
ALL_PROB_FLAG = true;

for i = 1:size(cst,1)

    for j = 1:numel(cst{i,6})
        if strcmp(cst{i,6}{j}.robustness,'PROB')
            PROB_FLAG = true;
        end

        if ~(strcmp(cst{i,6}{j}.robustness,'none') || strcmp(cst{i,6}{j}.robustness,'PROB'))
            ROB_FLAG = true;
        end

        if ~strcmp(cst{i,6}{j}.robustness, 'PROB')
            ALL_PROB_FLAG = false;
        end
    end
end
voiForOmegaIx = [];
%Set the structures to be included for prob calculation
if PROB_FLAG
    for i = 1:size(cst,1)
        for j=1:size(cst{i,6},2)
            if isa(cst{i,6}{j},'OmegaObjectives.matRad_TotalVariance')
                voiForOmegaIx = [voiForOmegaIx i];
            end
        end
    end
    voiForOmegaIx = unique(voiForOmegaIx);
end

%set scenarios to be included
if ROB_FLAG
    if ~isempty(allScen)
        %Get linear map of the scenario indexes
        [linMap(:,1), linMap(:,2), linMap(:,3)] = ind2sub([pln.multScen.numOfCtScen, pln.multScen.totNumShiftScen, pln.multScen.totNumRangeScen], allScen);

        %Select those that are not excluded by scen4D
        ixSelected4D = ismember(linMap(:,1),scen4D);

        %Select the linear indexes on dij correspoinding to those scenarios
        useScen = sub2ind([pln.multScen.numOfCtScen, pln.multScen.totNumShiftScen, pln.multScen.totNumRangeScen], linMap(ixSelected4D,1),linMap(ixSelected4D,2),linMap(ixSelected4D,3));

    else
        matRad_cfg.dispError('Trying to set robustness different from PROB but no scenarios have been stored');
    end
end

% If there are no scenarios and all structures have PROB robustness and no
% other kind of robustness, disable scenario calculation in optimization
% problem. This avoids the call to projectSingleScenario in BP
if (isempty(allScen) && PROB_FLAG && ~ROB_FLAG) || ALL_PROB_FLAG
    useScen = [];
end

% If dij scenarios are provided, the dose distribution will be computed for
% all the nominal CT scenarios by default.

if PROB_FLAG && ~((isfield(dij, 'physicalDoseExp') &&  isfield(dij, 'physicalDoseOmega')) || isfield(dij, 'mAlphaDoseExp') &&  isfield(dij, 'mAlphaDoseOmega'))
    [dij] = matRad_calculateProbabilisticQuantities(dij,cst,pln);
end

% one robust objective is defined
% Old robust opt
% % % %Check how to use 4D data
% % % if isfield(pln,'propOpt') && isfield(pln.propOpt,'scen4D')
% % %     scen4D = pln.propOpt.scen4D;
% % % else
% % %     scen4D = 1; %Use only first 4D scenario for optimization
% % % end
% % % 
% % % %If "all" provided, use all scenarios
% % % if isequal(scen4D,'all')
% % %     scen4D = 1:size(dij.physicalDose,1);
% % % end
% % % 
% % % if ~isempty(scen4D)
% % %     linMap = ismember(pln.multScen.linearMask(:,1),scen4D);
% % %     linIxDIJ = sub2ind([pln.multScen.numOfCtScen, pln.multScen.totNumShiftScen, pln.multScen.totNumRangeScen],  pln.multScen.linearMask(find(linMap),1), pln.multScen.linearMask(find(linMap),2), pln.multScen.linearMask(find(linMap),3))';%find(~cellfun(@isempty,dij.physicalDose(scen4D,:,:)))';
% % %     %Only select the indexes of the nominal ct Scenarios
% % %     linIxDIJ_nominalCT = sub2ind([pln.multScen.numOfCtScen, pln.multScen.totNumShiftScen, pln.multScen.totNumRangeScen], scen4D', ones(numel(scen4D),1), ones(numel(scen4D),1))';%find(~cellfun(@isempty,dij.physicalDose(scen4D,1,1)))';
% % % else
% % %     linMap = [];
% % %     linIxDIJ = [];
% % %     linIxDIJ_nominalCT = [];
% % % end
% % % 
% % % FLAG_CALC_PROB = false;
% % % 
% % % FLAG_ROB_OPT   = false;
% % % 
% % % FLAG_PROB_OPT = false;
% % % 
% % % for i = 1:size(cst,1)
% % %     for j = 1:numel(cst{i,6})
% % %         if strcmp(cst{i,6}{j}.robustness,'PROB') && numel(linIxDIJ) > 1
% % %             if (isfield(dij, 'physicalDoseExp') &&  isfield(dij, 'physicalDoseOmega')) || isfield(dij, 'mAlphaDoseExp') &&  isfield(dij, 'mAlphaDoseOmega')
% % %                 FLAG_CALC_PROB = true;
% % %             else
% % %                 matRad_cfg.dispWarning('Probabilistic quantities required for optimization but not present in dij.');
% % % 
% % %             end
% % %         end
% % %         if ~strcmp(cst{i,6}{j}.robustness,'none') && numel(linIxDIJ) > 1
% % %             FLAG_ROB_OPT = true;
% % % 
% % %         end
% % % 
% % %         if strcmp(cst{i,6}{j}.robustness, 'PROB') && isfield(dij, 'physicalDoseExp')
% % %             FLAG_PROB_OPT = true;
% % %         end
% % %     end
% % % end
% % % 
% % % if FLAG_CALC_PROB && ~((isfield(dij, 'physicalDoseExp') &&  isfield(dij, 'physicalDoseOmega')) || isfield(dij, 'mAlphaDoseExp') &&  isfield(dij, 'mAlphaDoseOmega'))
% % %     [dij] = matRad_calculateProbabilisticQuantities(dij,cst,pln);
% % % end
% % % 
% % % 
% % % 
% % % %This has to be fixed
% % % if isempty(linIxDIJ) && FLAG_PROB_OPT
% % %     scen4D = [];
% % %     linIxDIJ_nominalCT = [1:numel(dij.physicalDoseExp)];
% % % end
% % % % set optimization options
% % % if ~FLAG_ROB_OPT || FLAG_CALC_PROB     % if multiple robust objectives are defined for one structure then remove FLAG_CALC_PROB from the if clause
% % %     ixForOpt = scen4D;
% % % else
% % %     ixForOpt = linIxDIJ;
% % % 
% % % end
% % % 
% % % 
% % % voiIx = [];
% % % for i = 1:size(cst,1)
% % %     for j=1:size(cst{i,6},2)
% % %         if isa(cst{i,6}{j},'OmegaObjectives.matRad_TotalVariance')
% % %             voiIx = [voiIx i];
% % %         end
% % %     end
% % % end
% % % voiIx = unique(voiIx);

switch pln.bioParam.quantityOpt
    case 'effect'
        backProjection = matRad_EffectProjection;
    case 'RBExD'
        %Capture special case of constant RBE
        if strcmp(pln.bioParam.model,'constRBE')
            backProjection = matRad_ConstantRBEProjection;
        else
            backProjection = matRad_VariableRBEProjection;
        end
    case 'physicalDose'
        backProjection = matRad_DoseProjection;
    otherwise
        warning(['Did not recognize bioloigcal setting ''' pln.probOpt.bioOptimization '''!\nUsing physical dose optimization!']);
        backProjection = matRad_DoseProjection;
end

backProjection.scenarios    = useScen;

%Need to filter out the probabilities. Get a 3D mask with the
%probabilities, then select onlz the ones in useScen
maskScenProb = pln.multScen.scenMask;

maskScenProb(find(maskScenProb)) = pln.multScen.scenWeight; %pln.multScen.scenProb;
matRad_cfg.dispWarning('!!! Using scen weights and not probabilities !');

backProjection.scenarioProb = maskScenProb(useScen);

backProjection.nominalCtScenarios = nominalCTScen;
backProjection.useStructsForOmega = voiForOmegaIx;

%Older code
%Give scenarios used for optimization
% backProjection.scenarios    = ixForOpt;
% backProjection.scenarioProb = pln.multScen.scenProb;%(ixForOpt);%./sum(pln.multScen.scenProb);%pln.multScen.scenProb;
% backProjection.nominalCtScenarios = linIxDIJ_nominalCT;
% 
% backProjection.useStructsForOmega = voiForOmegaIx;



if ~isfield(pln.propOpt, 'visualizeSingleObjectives')
   pln.propOpt.visualizeSingleObjectives = matRad_cfg.propOpt.defaultVisualizeSingleObjectives;
end

% for i=1:size(cst)
%     for j=1:numel(cst{i,6})
%         objective = cst{i,6}{j};
%         if isa(objective,'DoseObjectives.matRad_DoseObjective') || isa(objective,'OmegaObjectives.matRad_OmegaObjective')
%             if isempty(objective.isActive)
% 
%             end
%         end
%     end
% end

optiProb = matRad_OptimizationProblemVisualization(backProjection);

optiProb.instantiateVisualization(cst);

optiProb.graphicOutput.active = pln.propOpt.visualizeSingleObjectives;


optiProb.quantityOpt = pln.bioParam.quantityOpt;
if isfield(pln,'propOpt') && isfield(pln.propOpt,'useLogSumExpForRobOpt')
    optiProb.useLogSumExpForRobOpt = pln.propOpt.useLogSumExpForRobOpt;
end

%Get Bounds
if ~isfield(pln.propOpt,'boundMU')
    pln.propOpt.boundMU = false;
end

if pln.propOpt.boundMU
    if (isfield(dij,'minMU') || isfield(dij,'maxMU')) && ~isfield(dij,'numParticlesPerMU')
        matRad_cfg.dispWarning('Requested MU bounds but number of particles per MU not set! Bounds will not be enforced and standard [0,Inf] will be used instead!');
    elseif ~isfield(dij,'minMU') && ~isfield(dij,'maxMU')
        matRad_cfg.dispWarning('Requested MU bounds but machine bounds not defined in dij.minMU & dij.maxMU! Bounds will not be enforced and standard [0,Inf] will be used instead!');
    else
        if isfield(dij,'minMU')
            optiProb.minimumW = dij.numParticlesPerMU .* dij.minMU / 1e6;
            matRad_cfg.dispInfo('Using lower MU bounds provided in dij!\n')
        end

        if isfield(dij,'maxMU')
            optiProb.maximumW = dij.numParticlesPerMU .* dij.maxMU / 1e6;
            matRad_cfg.dispInfo('Using upper MU bounds provided in dij!\n')
        end
    end
else
    matRad_cfg.dispInfo('Using standard MU bounds of [0,Inf]!\n')
end


if ~isfield(pln.propOpt,'optimizer')
    pln.propOpt.optimizer = 'IPOPT';
end

switch pln.propOpt.optimizer
    case 'IPOPT'
        optimizer = matRad_OptimizerIPOPT;
    case 'fmincon'
        optimizer = matRad_OptimizerFmincon;
    otherwise
        warning(['Optimizer ''' pln.propOpt.optimizer ''' not known! Fallback to IPOPT!']);
        optimizer = matRad_OptimizerIPOPT;
end
        
if ~optimizer.IsAvailable()
    matRad_cfg.dispError(['Optimizer ''' pln.propOpt.optimizer ''' not available!']);
end

optimizer = optimizer.optimize(wInit,optiProb,dij,cst);

wOpt = optimizer.wResult;
info = optimizer.resultInfo;



resultGUI = matRad_calcCubes(wOpt,dij);
resultGUI.wUnsequenced = wOpt;
resultGUI.usedOptimizer = optimizer;
resultGUI.info = info;

resultGUI.info.timePerIteration = resultGUI.info.cpu/resultGUI.info.iter;

for i=1:numel(optiProb.graphicOutput.data.objectiveFunctions)
    functionName = optiProb.graphicOutput.leg{i};
    resultGUI.costFunctions(i).name = functionName;
    resultGUI.costFunctions(i).values = optiProb.graphicOutput.data.objectiveFunctions(i).values;

end
resultGUI.costFunctions(i+1).name = 'total Function';
resultGUI.costFunctions(i+1).values = optiProb.graphicOutput.data.totFValues;

%Robust quantities
%if FLAG_ROB_OPT || numel(ixForOpt) > 1
if ROB_FLAG

    Cnt = 1;
    for i = find(~cellfun(@isempty,dij.physicalDose))'
        tmpResultGUI = matRad_calcCubes(wOpt,dij,i);
        resultGUI.([pln.bioParam.quantityVis '_' num2str(Cnt,'%d')]) = tmpResultGUI.(pln.bioParam.quantityVis);
        Cnt = Cnt + 1;
    end
end



% unblock mex files
clear mex
