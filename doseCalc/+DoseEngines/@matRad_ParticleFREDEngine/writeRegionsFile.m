function writeRegionsFile(this,fName, stf)
    matRad_cfg = MatRad_Config.instance();
    fred_cfg = MatRad_FREDConfig.instance();
    fID = fopen(fName, 'w');
    

    try
        fprintf(fID,'region<\n');
        fprintf(fID,'\tID=Phantom\n');
        
        if this.useWaterPhantom
            % TODO: check these coordinates dimensions
            fprintf(fID, '\tL=[%2.2f,%2.2f,%2.2f]\n', this.doseGrid.dimensions(2)*this.doseGrid.resolution.x/10,...
                                                      this.doseGrid.dimensions(1)*this.doseGrid.resolution.y/10,...
                                                      this.doseGrid.dimensions(3)*this.doseGrid.resolution.z/10);
            fprintf(fID, '\tvoxels=[%i,%i,%i]\n',  this.doseGrid.dimensions(2),  this.doseGrid.dimensions(1),  this.doseGrid.dimensions(3));
            fprintf(fID, '\tmaterial=water78\n');
        else
            fprintf(fID,'\tCTscan=inp/regions/%s\n', fred_cfg.patientFilename);
        end
            fprintf(fID,'\tO=[%i,%i,%i]\n', 0,0,0);
            fprintf(fID,'\tpivot=[0.5,0.5,0.5]\n');

        %This is gantry angle = 0
        fprintf(fID, '\tl=[%1.1f,%1.1f,%1.1f]\n', 1,0,0);
        fprintf(fID, '\tu=[%1.1f,%1.1f,%1.1f]\n', 0,-1,0);

        %if numel(this.scorers)>1
        switch this.currentVersion
            case '3.69.14'
                 if this.calcDoseDirect
                    fprintf(fID,'\tscore=[');
                else
                    fprintf(fID,'\tscoreij=[');
                end

            otherwise
                fprintf(fID,'\tscore=[');
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
        
        if ~this.useWaterPhantom
            switch this.HUtable
                case 'internal'
                    fprintf(fID, 'lUseInternalHU2Mat=t\n');
                otherwise
                    fprintf(fID, 'include: inp/regions/hLut.inp\n');
                    writeHlut(this.HUtable);
            end
            
            if this.HUclamping
                fprintf(fID, 'lAllowHUClamping=t\n');
            end        
        
        end

        if this.useWaterPhantom
            fprintf(fID, 'material<\n');
            fprintf(fID, '\tID=water78\n');
            fprintf(fID, '\tbasedOn=water\n');
            fprintf(fID, '\trho=1\n');
            fprintf(fID, '\tIpot=78\n');
            fprintf(fID, 'material>\n');
        end

    catch
        matRad_cfg.dispError('Failed to write regions file');
    end

    fclose(fID);
end

function writeHlut(hLutFile)
    fileName = fullfile('MCrun/inp/regions/', 'hLut.inp');
    template = fileread(hLutFile);

    newLut = fopen(fileName, 'w');
    fprintf(newLut, template);
    fclose(newLut);
end
