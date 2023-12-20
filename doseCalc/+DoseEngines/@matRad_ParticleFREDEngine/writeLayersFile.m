function writeLayersFile(this, fName, stf)
    matRad_cfg = MatRad_Config.instance();

    fID = fopen(fName, 'w');
    
    % This is hard coded here, could find a cleaner solution

    try
        fprintf(fID, 'include: inp/plan/beamlets.inp\n');
        currLayerId = 0;
        beamletCounter = 0;
         for i=1:numel(stf)
            for j=1:numel(stf(i).energies)
                currLayerId = currLayerId+1;
                fprintf(fID, 'def: L%i = {', currLayerId);
                fprintf(fID, '''Energy'': %i,', stf(i).energies(j));
                fprintf(fID, '''FWHM'': %3.2f,', stf(i).FWHMs(j));

                fprintf(fID, ' ''beamlets'': [');
                for k=1:stf(i).energyLayer(j).nBixels-1
                    fprintf(fID, 'S%i,',beamletCounter+k);
                end
                fprintf(fID, 'S%i]',beamletCounter+stf(i).energyLayer(j).nBixels);
                beamletCounter = beamletCounter+stf(i).energyLayer(j).nBixels;

            fprintf(fID, '}\n');
            end
         end
      
    catch
        matRad_cfg.dispError('Failed to write layers file');
    end

    fclose(fID);
end