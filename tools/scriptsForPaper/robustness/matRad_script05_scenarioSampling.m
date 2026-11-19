%% Script 5: Worst-case vs. randomly sampled uncertainty scenarios
%
% Demonstrates matRad's scenario classes (matRad/scenarios/) directly,
% with no ct/cst, no dose calculation, and no optimization involved.
% Robust planning in matRad represents setup and range uncertainty as a
% zero-mean Gaussian: independent setup shifts in x/y/z (SD = shiftSD,
% [mm]) and independent absolute/relative range errors (SD = rangeAbsSD
% [mm] / rangeRelSD [%]). A matRad_ScenarioModel subclass turns that
% uncertainty model into a concrete, finite set of scenarios:
%   - matRad_WorstCaseScenarios: one shift at a time is pushed to its
%     +/-wcSigma*SD extreme (a small, deterministic grid used for
%     worst-case robust optimization)
%   - matRad_RandomScenarios: scenarios are randomly drawn from the same
%     Gaussian (used for probabilistic robust optimization / uncertainty
%     quantification)
%
% This script builds one of each and plots the underlying 2D (bivariate)
% probability density of the x/y setup shifts, together with the
% scenarios actually sampled from it.

%% set matRad runtime configuration
clear;
matRad_rc;
matRad_cfg = MatRad_Config.instance();

%% Build the scenario models
% Neither class needs a ct - with no argument they default to a single
% (static, non-4D) ct scenario. All uncertainty-model properties
% (shiftSD, rangeAbsSD, rangeRelSD, wcSigma) are left at their matRad
% defaults and are shared between both models below, so the two are
% directly comparable.
scenWC = matRad_WorstCaseScenarios();
scenWC.wcSigma = 2;
scenWC.combinations = 'shift';

scenRand = matRad_RandomScenarios();
scenRand.nSamples = 30; % matRad's own default (10) is realistic for an actual robust optimization run, but too sparse to visualize a distribution - increased here for illustration only

% matRad_cfg.dispInfo('\n--- %s (%s) ---\n',scenWC.name,scenWC.shortName);
% scenWC.listAllScenarios();
% matRad_cfg.dispInfo('Total scenarios: %d\n',scenWC.totNumScen);
% 
% matRad_cfg.dispInfo('\n--- %s (%s) ---\n',scenRand.name,scenRand.shortName);
% matRad_cfg.dispInfo('Total scenarios: %d\n',scenRand.totNumScen);

%% Build the underlying 2D probability density (x/y setup shift)
% scenForProb columns are [ctScenIdx, xShift, yShift, zShift, absRangeShift, relRangeShift].
% x and y shifts are modeled as independent zero-mean Gaussians with SD =
% shiftSD(1) and shiftSD(2), i.e. a bivariate Gaussian with diagonal
% covariance. Computed by hand to avoid a Statistics Toolbox (mvnpdf)
% dependency.
sigmaX = scenWC.shiftSD(1); % [mm]
sigmaY = scenWC.shiftSD(2); % [mm]

xGrid = linspace(-4*sigmaX,4*sigmaX,200);
yGrid = linspace(-4*sigmaY,4*sigmaY,200);
[X,Y] = meshgrid(xGrid,yGrid);

% Mahalanobis distance
R = sqrt((X/sigmaX).^2 + (Y/sigmaY).^2);
pdf2D = exp(-0.5*R.^2) / (2*pi*sigmaX*sigmaY);

% 1/2/3 sigma ellipses (circles for sigmaX == sigmaY), drawn as plain lines
% so they can be exported with matRad_linesToTikzFile
sigmaLevels = [1 2 3];
phi = linspace(0,2*pi,200);

%% Plot the 2D density with sampled scenarios overlaid and export to tikz
% One figure per scenario model. Per model, saveDir/<shortName>/ gets:
%   - scenarioSampling.tex (+ png): the density colormap and colorbar
%     (exported via matlab2tikz, with the circles and dots hidden)
%   - contours.txt: the 1/2/3 sigma circles as \addplot commands
%   - scenarios.txt: the scenarios as an \addplot [only marks] list of dots
% Paste contours.txt and scenarios.txt inside the axis environment of the
% .tex file, after the \addplot graphics line.
saveDir = fullfile(matRad_cfg.primaryUserFolder, 'dataPaper', 'Figures', 'robustness', 'scenarioSampling');

savePlots = false;

modelsToPlot = {scenWC,'Worst case (wcScen)'; scenRand,'Random sampling (rndScen)'};
for i = 1:size(modelsToPlot,1)
    model = modelsToPlot{i,1};

    f = figure('Name',modelsToPlot{i,2},'Color','w');
    imagesc(xGrid,yGrid,pdf2D);
    axis xy equal tight;
    hold on;

    hCircles = gobjects(numel(sigmaLevels),1);
    for k = 1:numel(sigmaLevels)
        hCircles(k) = plot(sigmaLevels(k)*sigmaX*cos(phi),sigmaLevels(k)*sigmaY*sin(phi), ...
            '--','Color',[0.15 0.15 0.15],'LineWidth',1);
    end

    % scenarios projected onto the x/y plane - scenarios that only differ in
    % z shift or range error (e.g. the worst-case z/range scenarios, all at
    % x = y = 0) will stack on top of each other.
    scenX = model.scenForProb(:,2);
    scenY = model.scenForProb(:,3);
    hScen = scatter(scenX,scenY,50,'w','filled','MarkerEdgeColor','k','LineWidth',0.5);

    xlabel('x setup shift [mm]');
    ylabel('y setup shift [mm]');
    title(sprintf('%s: %d scenarios',modelsToPlot{i,2},model.totNumScen));
    cb = colorbar;
    cb.Label.String = 'probability density [mm^{-2}]';

    if savePlots
        figSaveDir = fullfile(saveDir, model.shortName);
        if ~exist(figSaveDir, 'dir')
            mkdir(figSaveDir);
        end

        % matRad_linesToTikzFile appends, so start from clean files
        contoursFile  = fullfile(figSaveDir, 'contours.txt');
        scenariosFile = fullfile(figSaveDir, 'scenarios.txt');
        if exist(contoursFile, 'file'),  delete(contoursFile);  end
        if exist(scenariosFile, 'file'), delete(scenariosFile); end

        % export only the density + colorbar through matlab2tikz, the
        % circles are written to contours.txt by matRad_saveFigToTikz
        title('');
        set([hCircles; hScen], 'Visible', 'off');
        matRad_saveFigToTikz(f, figSaveDir, 'scenarioSampling.tex', {hCircles});
        set([hCircles; hScen], 'Visible', 'on');

        % list of dots
        fid = fopen(scenariosFile, 'w');
        fprintf(fid, '%% %s: %d scenarios [xShift, yShift] in mm\n', model.name, model.totNumScen);
        fprintf(fid, '\\addplot [only marks, mark=*, mark size=\\scenarioMarkSize,fill=white, draw=black, forget plot] coordinates {');
        fprintf(fid, '(%.3f,%.3f) ', [scenX(:).'; scenY(:).']);
        fprintf(fid, '};\n');
        fclose(fid);
    end
end
