function [energySpread,E, costFunctionValue, exitFlag] = getEnergySpread(dicomFileName,machine,eIdx, nPeaks, initialGuess,verbose)
%% hard coded parametersmachine
    depthResolution = 0.1; % of topas simulation
    window = 40; % mm, size of window around peak to use for optimization
    meanEnergy = @(x) 11.39 * x^0.628 + 11.24;

    shifts = linspace(-3,3,nPeaks);%0.3*[-10:7]; % in mm, shift of the peeaks for optimization
    %initialGuess = 0.02;
%% Load dicom
    BP = dicomread(dicomFileName);
    BP = double(squeeze(BP));
    BP = smooth(squeeze(sum(BP, [2,3])));

%% Define resolution and depth vector of topas simulation
    x = [0:size(BP,1)-1]*depthResolution + depthResolution/2;

    
%% Cut central part

    [~,peakPosIdx] = max(BP);
    peakPos = x(peakPosIdx);
    if peakPos >window/1.5
        newDepthsIdx(1) = find(x<peakPos - window/1.5, 1, 'last');
    else
        newDepthsIdx(1) = 1;
    end
    newDepthsIdx(2) = find(x<peakPos + window/2, 1, 'last');
    cutDepths = x([newDepthsIdx(1):newDepthsIdx(2)]);
    cutPDD    = BP([newDepthsIdx(1):newDepthsIdx(2)]);

    if verbose >1
        figure;
        plot(cutDepths, cutPDD, 'o-');

    end

%% Oversample

    newDepths = linspace(cutDepths(1), cutDepths(end), 10*size(cutDepths,2));
    newPDD    = interp1(cutDepths, cutPDD, newDepths, 'spline');
    if verbose >1
        hold on;
        plot(newDepths, newPDD, '.-');
    end
%% Get range
    [maxV, maxI] = max(newPDD);

    [~, r80ind] = min(abs(newPDD(maxI:end) - 0.8 * maxV));
    r80ind = r80ind - 1;
    r80 = interp1(newPDD(maxI + r80ind - 1:maxI + r80ind + 1),newDepths(maxI + r80ind - 1:maxI + r80ind + 1), 0.8 * maxV);
    E = meanEnergy(r80);

%% Generate shifts
    shiftedR80 = r80 + shifts; %mm

    Ebins = arrayfun(@(r) meanEnergy(r), shiftedR80); % these are topas energies
    Erange = (Ebins(end)-Ebins(1))*12;
%% Get machine data

if verbose >1
    figure;
    for k =1:size(shiftedR80,2)
    
        plot(newDepths + shifts(k), newPDD./max(newPDD), '.-', 'color', 'k');
    
        hold on;
    end
    plot([machine.data(eIdx).depths], machine.data(eIdx).Z./max(machine.data(eIdx).Z), '.-', 'color', 'r');
     hold on;
     plot(newDepths, 0.5*ones(size(newDepths)), '-o');
end

%% Gen data for minimization
    optIntervall = [500:size(newDepths,2)-500];%[300:size(newDepths,2)-700];
    optDepths = newDepths(optIntervall);
    tPDDs = zeros(size(shifts,2),size(newDepths,2));
    for k=1:size(shifts,2)

        tPDDs(k,:) = interp1(newDepths+shifts(k), newPDD./max(newPDD), newDepths);
        tPDDs(k, isnan(tPDDs(k,:))) = 0;
    end
    PDDs = tPDDs(:,optIntervall);
    calcLinComb = @(w,bPs) w*bPs;

    sigmaMeV = initialGuess.*(meanEnergy(r80))/100;
    w0 = (1./(2*pi*sigmaMeV)).*exp((-(Ebins-meanEnergy(r80)).^2)./(2.*sigmaMeV.^2)); %ones(1,size(shifts,2));
    %w0 = exp((-(Ebins-meanEnergy(r80)).^2)./(2.*sigmaMeV.^2));
    w0 = w0/sum(w0);
    
    objective = interp1(machine.data(eIdx).depths, machine.data(eIdx).Z,newDepths,'spline');
    objective = objective(optIntervall)./max(objective);
    %objective = newPDD(optIntervall)./max(newPDD);
    if verbose >0
        figure;
        plot(optDepths, objective, '--');
        hold on;
        plot(optDepths, calcLinComb(w0,PDDs));
%         for k = size(shifts,2)
%             plot(optDepths, PDDs(k,:),'.-');
%             plot(optDepths, w0(k).*PDDs(k,:), '--');
%         end
        figure;
        plot(Ebins, w0, '.-');
    end
    
    
    options = optimoptions('lsqcurvefit');
    numberfOfVariables = size(w0,2);
    
    MaxEvaluatuions = 200;
    options.MaxFunctionEvaluations = MaxEvaluatuions*numberfOfVariables;

    [resW, costFunctionValue, ~, exitFlag] = lsqcurvefit(@(w,D) calcLinComb(w,D), w0, PDDs,objective, zeros(size(w0)), ones(size(w0)), options);

    resultPDD = calcLinComb(resW,PDDs);
    if verbose >0
        figure;
        plot(optDepths,resultPDD,'.-');
        hold on;
        plot(optDepths, objective, '--');
        for k = 1:size(shifts,2)
            plot(optDepths, resW(k).*PDDs(k,:), '--');
        end
         plot(machine.data(eIdx).depths,machine.data(eIdx).Z./max(machine.data(eIdx).Z), '.-');
         xlim([0,r80+0.3*r80]);
    end

%% Plot weights

    if verbose >0
        figure;
        plot(Ebins, resW./sum(resW), '-o');
    end
%% Gauss fit on weights
    resF = fit(Ebins', resW', 'gauss1');
    sigma = resF.c1*sqrt(1/2);
    energySpread = (sigma*100)/(resF.b1);
%% Plot fit
if verbose >0
    energiesForPlot = linspace(Ebins(1)-10, Ebins(end)+10, 100*size(Ebins,2));
    hold on;
    plot(energiesForPlot, resF(energiesForPlot), '.-');

end



%% Replot
% newResultPDD = calcLinComb(resF(Ebins)',PDDs);
% if verbose >0
%         figure;
%         plot(optDepths,newResultPDD,'.-');
%         hold on;
%         plot(optDepths, objective, '--');
%         for k = 1:size(shifts,2)
%             plot(optDepths, resF(Ebins(k)).*PDDs(k,:), '--');
%         end
%          plot(machine.data(eIdx).depths,machine.data(eIdx).Z./max(machine.data(eIdx).Z), '.-');
%          xlim([0,r80+0.3*r80]);
% end
end