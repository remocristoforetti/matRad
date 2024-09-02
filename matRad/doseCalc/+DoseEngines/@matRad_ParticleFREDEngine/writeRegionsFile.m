function writeRegionsFile(this,fName, stf)
    matRad_cfg = MatRad_Config.instance();
    fID = fopen(fName, 'w');
    
    try
        
        fprintf(fID,'region<\n');
        fprintf(fID,'\tID=Phantom\n');
        fprintf(fID,'\tCTscan=inp/regions/%s\n', this.patientFilename);
        fprintf(fID,'\tO=[%i,%i,%i]\n', 0,0,0);
        fprintf(fID,'\tpivot=[0.5,0.5,0.5]\n');

        %This is gantry angle = 0
        fprintf(fID, '\tl=[%1.1f,%1.1f,%1.1f]\n', 1,0,0);
        fprintf(fID, '\tu=[%1.1f,%1.1f,%1.1f]\n', 0,-1,0);

        if this.calcDoseDirect
            fprintf(fID,'\tscore=[');
        else
            fprintf(fID,'\tscoreij=[');
        end

        if numel(this.scorers)>1

            for k=1:size(this.scorers,2)-1
                fprintf(fID,'%s,', this.scorers{k});
            end
        end
        fprintf(fID,'%s]\n', this.scorers{end});    
        
        fprintf(fID,'region>\n');

        fprintf(fID, 'region<\n');
        fprintf(fID, '\tID=Room\n');
        fprintf(fID, '\tmaterial=%s\n', this.roomMaterial);
        fprintf(fID, 'region>\n');
        
        switch this.HUtable
            case 'internal'
                fprintf(fID, 'lUseInternalHU2Mat=t\n');
            otherwise
                fprintf(fID, 'include: inp/regions/hLut.inp\n');
                this.writeHlut(this.HUtable);
        end
        
        if this.HUclamping
            fprintf(fID, 'lAllowHUClamping=t\n');
        end

    catch ME
        matRad_cfg.dispError(['Failed to write regions file. Exit with error: ', ME.message]);
    
    end

    fclose(fID);
end


