function matRad_writeTickzDVH(dvh, filename)
    
    fid = fopen(filename, 'w');
    if fid == -1
        error('Could not open file "%s" for writing.', filename);
    end

    % Write PGFPlots header
    fprintf(fid, '\\addplot [color=mycolor1, line width=2.0pt, forget plot]');
    fprintf(fid, '  table[row sep=crcr]{%\n');

    % Write data rows
    for i = 1:length(dvh.doseGrid)
        fprintf(fid, '%.6f\t%.6f\\\\ ', dvh.doseGrid(i), dvh.volumePoints(i));
    end

    % Close the table and file
    fprintf(fid, '};\n');
    fclose(fid);

    fprintf('Data successfully written to "%s".\n', filename);

end