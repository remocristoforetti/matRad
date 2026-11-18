function matRad_writeTickzDVH_Fill(dvh_h,dvh_l, filename)
    
    fid = fopen(filename, 'w');
    if fid == -1
        error('Could not open file "%s" for writing.', filename);
    end

    % Write PGFPlots header
    fprintf(fid, '\\addplot [area legend, draw=none, fill=mycolor2, fill opacity=0.2, forget plot]');
    fprintf(fid, '  table[row sep=crcr]{%\n');

    % Write data rows
    for i = 1:length(dvh_h.doseGrid)
        fprintf(fid, '%.6f\t%.6f\\\\ ', dvh_h.doseGrid(i), dvh_h.volumePoints(i));
    end
    
    f_doseGrid = flip(dvh_l.doseGrid);

    f_volume   = flip(dvh_l.volumePoints);

    for i = 1:length(dvh_l.doseGrid)
        fprintf(fid, '%.6f\t%.6f\\\\ ', f_doseGrid(i),f_volume(i));
    end

    % Close the table and file
    fprintf(fid, '};\n');
    fclose(fid);

    fprintf('Data successfully written to "%s".\n', filename);

end