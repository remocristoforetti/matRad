function matRad_linesToTikzFile(hLines, filename, opts, lineWidthName)
% LINESTOTIKZFILE Write MATLAB Line objects as pgfplots \addplot commands.
%
%   linesToTikzFile(hContour, filename)
%   linesToTikzFile(hContour, filename, opts)
%
%   hLines   - cell array of Line objects (e.g. hContour from your custom
%              contour-drawing function), an array of Line objects, or a
%              single Line object. Color, LineStyle, LineWidth and Marker
%              are read straight off each object.
%   filename - output .txt path
%   opts     - optional struct:
%       .forgetPlot     logical, add 'forget plot' so lines are skipped in
%                        the legend (default true)
%       .precision      decimal digits for coordinates (default 3)
%       .colorPrecision decimal digits for rgb color components (default 4)
%       .namePrefix     prefix used for \definecolor names (default 'lineColor')
%       .append         logical, append to filename instead of overwriting
%                        it, so repeated calls accumulate in the same file
%                        (default true)
%
%   Paste the resulting file's contents inside the tikz axis environment,
%   after the \addplot graphics [...] {...}; line(s).

if nargin < 3, opts = struct(); end
if ~isfield(opts,'forgetPlot'),     opts.forgetPlot = true; end
if ~isfield(opts,'precision'),      opts.precision = 3; end
if ~isfield(opts,'colorPrecision'), opts.colorPrecision = 4; end
if ~isfield(opts,'namePrefix'),     opts.namePrefix = 'lineColor'; end
if ~isfield(opts,'append'),         opts.append = true; end

if ~exist('lineWidthName', 'var') || isempty(lineWidthName)
    lineWidthName = 'lineWidth';
end

if iscell(hLines)
    lines = hLines;
else
    lines = num2cell(hLines);
end

if ~isfield(opts,'lineWidth')
    opts.lineWidth = lines{1}.LineWidth;
end

if opts.append
    fid = fopen(filename,'a');
else
    fid = fopen(filename,'w');
end
if fid == -1
    error('linesToTikzFile:fopen', 'Could not open %s for writing.', filename);
end
cleanupObj = onCleanup(@() fclose(fid));

% fprintf(fid, '\\def\\%s{%.2fpt}\n\n', lineWidthName, opts.lineWidth);

coordFmt = sprintf('%%.%df', opts.precision);
colorFmt = sprintf('%%.%df', opts.colorPrecision);

styleMap = containers.Map( ...
    {'-', '--', ':', '-.', 'none'}, ...
    {'solid', 'dashed', 'dotted', 'dashdotted', 'solid'});

markerMap = containers.Map( ...
    {'none','o','+','*','x','square','diamond','^','v','>','<','pentagram','hexagram','.'}, ...
    {'none','o','+','asterisk','x','square','diamond','triangle','triangle','triangle','triangle','star','star','*'});

for i = 1:numel(lines)
    ln = lines{i};
    x = ln.XData;
    y = ln.YData;
    if isempty(x)
        continue;
    end

    colorName = sprintf('%s%03d%03d%03d', opts.namePrefix, ...
        round(ln.Color(1)*255), round(ln.Color(2)*255), round(ln.Color(3)*255));
    fprintf(fid, ['\\definecolor{%s}{rgb}{' colorFmt ',' colorFmt ',' colorFmt '}\n'], ...
        colorName, ln.Color(1), ln.Color(2), ln.Color(3));

    if isKey(styleMap, ln.LineStyle)
        tikzStyle = styleMap(ln.LineStyle);
    else
        tikzStyle = 'solid';
    end

    if isKey(markerMap, ln.Marker)
        tikzMarker = markerMap(ln.Marker);
    else
        tikzMarker = 'none';
    end

    plotOpts = sprintf('color=%s, %s, line width=\\%s, mark=%s', ...
        colorName, tikzStyle,lineWidthName, tikzMarker);
    if opts.forgetPlot
        plotOpts = [plotOpts, ', forget plot']; %#ok<AGROW>
    end

    fprintf(fid, '%% line %d: %d points, color=[%.3f %.3f %.3f]\n', ...
        i, numel(x), ln.Color(1), ln.Color(2), ln.Color(3));
    fprintf(fid, '\\addplot [%s] coordinates {', plotOpts);
    fprintf(fid, ['(' coordFmt ',' coordFmt ') '], [x(:).'; y(:).']);
    fprintf(fid, '};\n\n');
end

end
