function writePlanFile(this,fName, stf)
    matRad_cfg = MatRad_Config.instance();

    fID = fopen(fName, 'w');
    
    try
        fprintf(fID, 'include: inp/plan/fields.inp\n');
        fprintf(fID, 'def: plan = {''SAD'': %2.3f, ''Fields'': [', stf(1).BAMStoIsoDist);
        for i=1:numel(stf)-1
            fprintf(fID,'F%i, ', i-1);
        end
        fprintf(fID, 'F%i]}', numel(stf)-1);
    catch
        matRad_cfg.dispError('Failed to write plan file');
    end

    fclose(fID);

end