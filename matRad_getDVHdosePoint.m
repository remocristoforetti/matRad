function [dPoint,dStd] = matRad_getDVHdosePoint(dvhs, structIdx, D_point, doseGrid, weights, fractions, upperLimit, lowerLimit)

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
        upLim = D_point+upperLimit;
        
        lowLim = D_point-lowerLimit;
        
        tmpDVHLimits(1) = find(currVolumePoints < upLim, 1,'first')-1;
        tmpDVHLimits(2) = find(currVolumePoints < lowLim, 1,'first');
    
        [currVolumePointsInterp, currVidx] = unique(currVolumePoints(tmpDVHLimits(1):tmpDVHLimits(2)));
        currDosePoints = doseGrid(tmpDVHLimits(1):tmpDVHLimits(2));
        currDosePoints = currDosePoints(currVidx);
        currPoints(scenIdx) = interp1(currVolumePointsInterp, currDosePoints*fractions,D_point);
    end

    dPoint = sum(weights'.*currPoints);
    dStd = std(currPoints,weights');
end