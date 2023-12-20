function writeFieldsFile(this,fName,stf)

    matRad_cfg = MatRad_Config.instance();

    fID = fopen(fName, 'w');
    
    % This is hard coded here, could find a cleaner solution


    %Notes
    % The f vector in field is the direction of propagation (by
    % FRED)
    % It is referred to the e3 element of the base, this is for the
    % specific component, thus L(3) is the size of the field in the
    % direction fo propagation. O is referred to the room
    % coordinates system
    layerCounter = 0;
    try
        fprintf(fID, 'include: inp/plan/layers.inp\n');
        fprintf(fID, 'def: nFields = %i\n',numel(stf));
        for i=1:numel(stf)
            currFieldID      = i-1;
            fieldLimX = max(abs([stf(i).energyLayer.rayPosX])) + 10*max([stf(i).FWHMs]);
            fieldLimY = max(abs([stf(i).energyLayer.rayPosY])) + 10*max([stf(i).FWHMs]);

            fprintf(fID, 'def: F%i = {', currFieldID);
            fprintf(fID, '''fieldNumber'': %i,', currFieldID);
            fprintf(fID, '''GA'': %3.2f,', stf(i).gantryAngle);
            fprintf(fID, ' ''CA'': %3.2f,', stf(i).couchAngle);

            fprintf(fID, ' ''ISO'': [%3.2f,%3.2f,%3.2f], ',stf(i).isoCenter);
            fprintf(fID, ' ''dim'': [%3.2f,%3.2f,%3.2f], ',[fieldLimX, fieldLimY,0.1]);

            fprintf(fID, ' ''Layers'': [');
            for j=1:numel(stf(i).energies)-1
                fprintf(fID, 'L%i,',layerCounter+j);
            end
            fprintf(fID, 'L%i]',layerCounter+numel(stf(i).energies));
            layerCounter = layerCounter+numel(stf(i).energies);

            fprintf(fID, '}\n');

        end

    catch
        matRad_cfg.dispError('Failed to write field file');
    end

    fclose(fID);
end