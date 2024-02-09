function writeRunFile(~, fName)
    matRad_cfg = MatRad_Config.instance();

    fID = fopen(fName, 'w');

    try
        %fprintf(fID, 'include: inp/funcs.inp\n');
        fprintf(fID, 'include: inp/regions/regions.inp\n');
        fprintf(fID, 'include: inp/plan/planDelivery.inp\n');
    catch
        matRad_cfg.dispError('Failed to write run file');
    end

    fclose(fID);
end