function [dPoint,dStd,allPoints, weights] = matRad_getDVHVolumePoint(dvhs, structIdx, V_point, doseGrid, weights, fractions, upperLimit, lowerLimit)

    matRad_cfg = MatRad_Config.instance();
    if ~exist('doseGrid', 'var')
        matRad_cfg.dispError('No dose grid defined');
    end

     if ~exist('fractions', 'var') || isempty(fractions)
        fractions = 1;
     end

    if ~exist('upperLimit', 'var') || isempty(upperLimit)

        upperLimit = 0;%15;
    end

    if ~exist('lowerLimit', 'var') || isempty(lowerLimit)
        lowerLimit = 0;%3;
    end

    volemePointsData = dvhs(structIdx).volumePoints;
    totNumScen = size(volemePointsData,1);

    if ~exist('weights', 'var') || isempty(weights)
        weights = ones(totNumScen,1)./totNumScen;
    end

    if ~eq(round(sum(weights)*10^10)/10^10 , 1)
        matRad_cfg.dispWarning('Weights are not normalized. Normalizing.');
        weights = weights./sum(weights);
    end

    for scenIdx=1:totNumScen
        currDVH = dvhs(structIdx);
        currVolumePoints = currDVH.volumePoints(scenIdx,:);
        upLim = V_point+upperLimit;
        
        lowLim = V_point-lowerLimit;
        
        tmpDVHLimits(2) = find(doseGrid > upLim, 1,'first')+1;
        tmpDVHLimits(1) = find(doseGrid < lowLim, 1,'last')-1;
    
        [currDosePointsInterp] = unique(doseGrid(tmpDVHLimits(1):tmpDVHLimits(2)));
        currVPoints = currVolumePoints(tmpDVHLimits(1):tmpDVHLimits(2));
        %currVPoints = currVolumePoints(currDidx);
        currPoints(scenIdx) = interp1(currDosePointsInterp, currVPoints,V_point);
    end

    dPoint = sum(weights'.*currPoints);
    dStd = std(currPoints,weights');
    allPoints = currPoints;
end