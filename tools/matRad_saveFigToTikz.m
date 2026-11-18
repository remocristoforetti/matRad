function matRad_saveFigToTikz(figureHandle, saveDir, fileName, hContours, hDoseLines)

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

matlab2tikz('figurehandle', figureHandle, 'width', '\textwidth', 'floatFormat', '%.3g', 'fileName', fullfile(saveDir, fileName));

if exist('hContours', 'var') && ~isempty(hContours)
    for i=1:numel(hContours)
        contoursFilename = fullfile(saveDir, 'contours.txt');
        matRad_linesToTikzFile(hContours{i}, contoursFilename, [],'contoursLineWidth');
    end
end


if exist('hDoseLines', 'var') && ~isempty(hDoseLines)

    for i=1:numel(hDoseLines)
        doseLinesFilename = fullfile(saveDir, 'doseLines.txt');
        matRad_linesToTikzFile(hDoseLines{i}, doseLinesFilename,[], 'isoDoseLineWidth');
    end
end

end