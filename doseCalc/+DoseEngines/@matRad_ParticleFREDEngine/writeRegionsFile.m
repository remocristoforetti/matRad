function writeRegionsFile(this,fName, stf)
    matRad_cfg = MatRad_Config.instance();
    fred_cfg = MatRad_FREDConfig.instance();
    fID = fopen(fName, 'w');
    
    try
        fprintf(fID,'region<\n');
        fprintf(fID,'\tID=Phantom\n');
        fprintf(fID,'\tCTscan=inp/regions/%s\n', fred_cfg.patientFilename);
        fprintf(fID,'\tO=[%i,%i,%i]\n', 0,0,0);
        fprintf(fID,'\tpivot=[0.5,0.5,0.5]\n');

        %This is gantry angle = 0
        fprintf(fID, '\tl=[%1.1f,%1.1f,%1.1f]\n', 1,0,0);
        fprintf(fID, '\tu=[%1.1f,%1.1f,%1.1f]\n', 0,-1,0);

        if numel(this.scorers)>1
            fprintf(fID,'\tscore=[');
            for k=1:size(this.scorers,2)-1
                fprintf(fID,'%s,', this.scorers{k});
            end
            fprintf(fID,'%s]\n', this.scorers{end});
        else
            fprintf(fID,'\tscore=[%s]\n', this.scorers{1});
        end

        fprintf(fID,'region>\n');

        fprintf(fID, 'region<\n');
        fprintf(fID, '\tID=Room\n');
        fprintf(fID, '\tmaterial=Air\n');
        fprintf(fID, 'region>\n');
        
        if this.useInternalHUConversion
            fprintf(fID, 'lUseInternalHU2Mat=t\n');
            if this.HUclamping

                fprintf(fID, 'lAllowHUClamping=t');
            end
        end
    catch
        matRad_cfg.dispError('Failed to write regions file');
    end

    fclose(fID);
end
