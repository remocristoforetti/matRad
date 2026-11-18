function matRad_writeTickzPlot(x,y, filename, counter)
    
    if ~exist('counter', 'var') || isempty(counter)
        counter = 1;
    end

    if counter>1
        fid = fopen(filename, 'a');
    else
        fid = fopen(filename, 'w');
    end

    if fid == -1
        error('Could not open file "%s" for writing.', filename);
    end

    % Write PGFPlots header
    fprintf(fid, '\\addplot [color=mycolor%d, line width=2.0pt, forget plot]', counter);
    fprintf(fid, '  table[row sep=crcr]{%\n');

    % Write data rows
    for i = 1:length(x)
        fprintf(fid, '%.6f\t%.6f\\\\ ', x(i), y(i));
    end

    % Close the table and file
    fprintf(fid, '};\n');
    fclose(fid);

    fprintf('Data successfully written to "%s".\n', filename);

end