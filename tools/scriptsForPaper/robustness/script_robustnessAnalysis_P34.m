%% Script: Robustness analysis of the nominal, STOCH and COWC plans (P34)
% 1) loads patient, plans and the nominal/robust dijs
% 2) builds the struct 'results' with one field per plan, holding the
%    nominal and all scenario dose distributions
% 3) analysis on 'results' (DVHs, DVH bands, QIs, voxel-wise stats, plots)
% Doses are total RBExD: physicalDose * dij.RBE * numOfFractions.
clear;
matRad_rc;
matRad_cfg = MatRad_Config.instance();

dijDir    = 'D:\r408i\P34wc2Scenarios';
planDir   = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','robustness','P34');
planNames = {'nominal','STOCH','COWC'};

saveResults = false; % save the 'results' struct with all distributions (large file)

%% Load patient, plans and dijs
load('Lung_Patient_34.mat','ct','cst'); % original cst, used for all evaluations
ct = matRad_getWorldAxes(ct);

nPlans = numel(planNames);
w = cell(1,nPlans);
for p = 1:nPlans
    tmp  = load(fullfile(planDir,[planNames{p} '.mat']),'wOpt','pln');
    w{p} = tmp.wOpt;
end
numOfFractions = tmp.pln.numOfFractions;

dijNom = load(fullfile(dijDir,'dij_nominal.mat'),'dij');
dijNom = dijNom.dij;

dijRob = load(fullfile(dijDir,'dij.mat'),'dij','pln');
multScen = dijRob.pln.multScen;
dijRob   = dijRob.dij;

doseGrid = dijNom.doseGrid;
cstGrid  = matRad_resizeCstToGrid(cst, ct.x,ct.y,ct.z, doseGrid.x,doseGrid.y,doseGrid.z);
nScen    = multScen.totNumScen;

%% Build results: nominal and scenario distributions per plan
calcDose = @(D,w) reshape(full(D*w)*dijNom.RBE*numOfFractions, doseGrid.dimensions);

results = struct();
for p = 1:nPlans
    matRad_cfg.dispInfo('Computing distributions for plan %s\n',planNames{p});

    res.w           = w{p};
    res.doseNominal = calcDose(dijNom.physicalDose{1},w{p});

    res.doseScen = cell(1,nScen);
    for s = 1:nScen
        lm = multScen.linearMask(s,:); % [ctScen shiftScen rangeScen]
        res.doseScen{s} = calcDose(dijRob.physicalDose{lm(1),lm(2),lm(3)},w{p});
    end

    results.(planNames{p}) = res;
end

% sanity check: nominal scenario of the robust dij vs. the nominal dij
diffNom = max(abs(results.(planNames{1}).doseScen{1}(:) - results.(planNames{1}).doseNominal(:)));
matRad_cfg.dispInfo('Max difference nominal dij vs. robust dij scenario 1 (%s plan): %.3g Gy(RBE)\n',planNames{1},diffNom);

clear dijNom dijRob;

% if saveResults
%     save(fullfile(planDir,'robustnessResults.mat'),'results','planNames','cstGrid','doseGrid','multScen','numOfFractions','-v7.3');
% end

%% Settings
nDvhPoints  = 1000;
plotSlice   = 63; % dose-grid slice for the nominal distribution plots
percentiles = [5,25,75,95];
qiNames     = {'Dmean','D2','D50','D95'};
structNames = cstGrid(:,2);
nStructs    = size(cstGrid,1);

plotStructs  = {'GTV','PTV','ESOPHAGUS','HEART','LUNGS_NOT_GTV','CORD'}; % structures shown in the DVH plot
saveAnalysis = true; % save the analysis (everything except the scenario dose cubes)

% one common DVH dose axis for all plans and scenarios
maxNominal  = max(cellfun(@(pn) max(results.(pn).doseNominal(:)), planNames));
dvhDoseGrid = linspace(0, 1.3*maxNominal, nDvhPoints);

%% DVHs, DVH bands, QIs and voxel-wise statistics
for p = 1:nPlans
    res = results.(planNames{p});

    % nominal
    res.dvhNominal = matRad_calcDVH(cstGrid, res.doseNominal, 'cum', dvhDoseGrid);

    % scenarios: DVHs [nStructs x nDvhPoints x nScen]
    res.dvhScen = zeros(nStructs,nDvhPoints,nScen);
    for s = 1:nScen
        dvh = matRad_calcDVH(cstGrid, res.doseScen{s}, 'cum', dvhDoseGrid);
        for i = 1:nStructs
            res.dvhScen(i,:,s) = dvh(i).volumePoints;
        end
    end

    % mean and percentile DVHs over scenarios, same struct layout as matRad_calcDVH
    for structIdx = 1:nStructs
        % [nScen x nDvhPoints]
        allScenariosVolumePoints = squeeze(res.dvhScen(structIdx,:,:)).';

        res.dvhMean(structIdx).doseGrid     = dvhDoseGrid;
        res.dvhMean(structIdx).volumePoints = mean(allScenariosVolumePoints,1);
        res.dvhMean(structIdx).name         = structNames{structIdx};

        for percentileP = percentiles
            res.(['dvhPercentile_', num2str(percentileP)])(structIdx).doseGrid     = dvhDoseGrid;
            res.(['dvhPercentile_', num2str(percentileP)])(structIdx).volumePoints = prctile(allScenariosVolumePoints,percentileP,1);
            res.(['dvhPercentile_', num2str(percentileP)])(structIdx).name         = structNames{structIdx};
        end
    end

    % DVH bands: envelope over all scenarios
    res.dvhBandLower = min(res.dvhScen,[],3);
    res.dvhBandUpper = max(res.dvhScen,[],3);

    % % QIs: res.qi.<qiName>.nominal [nStructs x 1], .scen [nStructs x nScen],
    % % .min/.max/.mean over scenarios [nStructs x 1]
    % qiNominal = zeros(nStructs,numel(qiNames));
    % qiScen    = zeros(nStructs,nScen,numel(qiNames));
    % for i = 1:nStructs
    %     qiNominal(i,:) = getDoseStats(res.doseNominal(cstGrid{i,4}{1}));
    %     for s = 1:nScen
    %         qiScen(i,s,:) = getDoseStats(res.doseScen{s}(cstGrid{i,4}{1}));
    %     end
    % end
    % for k = 1:numel(qiNames)
    %     res.qi.(qiNames{k}).nominal = qiNominal(:,k);
    %     res.qi.(qiNames{k}).scen    = qiScen(:,:,k);
    %     res.qi.(qiNames{k}).min     = min(qiScen(:,:,k),[],2);
    %     res.qi.(qiNames{k}).max     = max(qiScen(:,:,k),[],2);
    %     res.qi.(qiNames{k}).mean    = mean(qiScen(:,:,k),2);
    % end

    % voxel-wise statistics over scenarios
    allScen = cat(4,res.doseScen{:});
    res.doseMin  = min(allScen,[],4);
    res.doseMax  = max(allScen,[],4);
    res.doseMean = mean(allScen,4);

    results.(planNames{p}) = res;
    clear res allScen;
end

%% Save the analysis (results without the scenario dose cubes)
if saveAnalysis
    analysis = struct();
    for p = 1:nPlans
        analysis.(planNames{p}) = rmfield(results.(planNames{p}),'doseScen');
    end
    save(fullfile(planDir,'robustnessAnalysis.mat'),'analysis','planNames','structNames','qiNames', ...
        'percentiles','dvhDoseGrid','doseGrid','multScen','numOfFractions','-v7.3');
end

%% Plot nominal distributions (on the dose grid)
ctGridDose = ct;
ctGridDose.cubeDim    = doseGrid.dimensions;
ctGridDose.cubeHU     = {interp3(ct.x,ct.y',ct.z,ct.cubeHU{1},doseGrid.x,doseGrid.y',doseGrid.z)};
ctGridDose.cube       = ctGridDose.cubeHU;
ctGridDose.x          = doseGrid.x;
ctGridDose.y          = doseGrid.y;
ctGridDose.z          = doseGrid.z;
ctGridDose.resolution = doseGrid.resolution;
ctGridDose = rmfield(ctGridDose,intersect(fieldnames(ctGridDose),{'size','origin'}));

doseWindow = [0 maxNominal];
figure('Color','w','WindowState','maximized');
for p = 1:nPlans
    subplot(1,nPlans,p);
    matRad_plotSlice(ctGridDose,'axesHandle',gca,'cst',cstGrid,'dose',results.(planNames{p}).doseNominal, ...
        'plane',3,'slice',plotSlice,'doseWindow',doseWindow,'title',[planNames{p} ' - nominal']);
end

%% Plot nominal DVHs with percentile bands
% per structure: min-max envelope (lightest), 5-95 and 25-75 percentile
% bands, scenario mean (dashed) and nominal (solid)
plotStructIdx = find(ismember(structNames,plotStructs));
colors = lines(numel(plotStructIdx));
shade  = @(lower,upper,c,a) fill([dvhDoseGrid fliplr(dvhDoseGrid)], [lower fliplr(upper)], ...
    c, 'FaceAlpha', a, 'EdgeColor', 'none', 'HandleVisibility', 'off');

figure('Color','w','WindowState','maximized');
for p = 1:nPlans
    res = results.(planNames{p});
    subplot(1,nPlans,p);
    hold on;
    for j = 1:numel(plotStructIdx)
        i = plotStructIdx(j);
        shade(res.dvhBandLower(i,:), res.dvhBandUpper(i,:), colors(j,:), 0.10);
        shade(res.dvhPercentile_5(i).volumePoints,  res.dvhPercentile_95(i).volumePoints, colors(j,:), 0.20);
        shade(res.dvhPercentile_25(i).volumePoints, res.dvhPercentile_75(i).volumePoints, colors(j,:), 0.35);
        plot(dvhDoseGrid, res.dvhMean(i).volumePoints, '--', 'Color', colors(j,:), 'LineWidth', 1, ...
            'HandleVisibility', 'off');
        plot(dvhDoseGrid, res.dvhNominal(i).volumePoints, '-', 'Color', colors(j,:), 'LineWidth', 1.5, ...
            'DisplayName', strrep(structNames{i},'_',' '));
    end
    grid on; box on;
    xlim([0 dvhDoseGrid(end)]); ylim([0 100]);
    xlabel('RBExD [Gy(RBE)]');
    ylabel('volume [%]');
    title(sprintf('%s: nominal DVH + bands (%d scenarios)',planNames{p},nScen));
end
legend('show','Location','northeastoutside');

%% Local functions
% function stats = getDoseStats(doseInVoi)
% % [Dmean, D2, D50, D95] of the dose values in a structure
% if isempty(doseInVoi)
%     stats = NaN(1,4);
%     return;
% end
% sorted = sort(doseInVoi,'descend');
% n  = numel(sorted);
% Dx = @(x) sorted(max(1,ceil(x/100*n)));
% stats = [mean(doseInVoi), Dx(2), Dx(50), Dx(95)];
% end
