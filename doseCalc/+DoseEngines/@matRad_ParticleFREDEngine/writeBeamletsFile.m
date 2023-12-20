function writeBeamletsFile(this, fName, stf)
    matRad_cfg = MatRad_Config.instance();

    fID = fopen(fName, 'w');
    
    try
        fprintf(fID, 'nprim = %i\n', this.numHistoriesDirect);
        beamletCounter = 0;
        for i=1:numel(stf)
            for j=1:numel(stf(i).energies)
                for k=1:stf(i).energyLayer(j).nBixels
                    beamletCounter = beamletCounter+1;
                    fprintf(fID, 'def: S%i = {', beamletCounter);
                    fprintf(fID, '''beamletID'': %i,', beamletCounter);
                    fprintf(fID, '''P'': [%2.3f,%2.3f,%2.3f], ', stf(i).energyLayer(j).rayPosX(k),stf(i).energyLayer(j).rayPosY(k),0);
                    fprintf(fID, '''v'': [%2.5f,%2.5f,%i], ', stf(i).energyLayer(j).rayDivX(k),stf(i).energyLayer(j).rayDivY(k),1);
                    fprintf(fID, '''w'': %2.4f}\n',stf(i).energyLayer(j).numOfPrimaries(k));                            
                end
            end
        end
    catch
        matRad_cfg.dispError('Failed to write beamlets file');
    end

    fclose(fID);

end