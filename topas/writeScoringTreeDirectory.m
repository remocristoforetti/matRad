function writeScoringTreeDirectory(phantoms, Ions, topas_cfg, EParam, includePDD, includeDS, includeSP,includeEdEventSpectrum,includeProtonLET)

   wDir = topas_cfg.workingDir;
   %make Scorer directory
   if ~exist([wDir,filesep,'Scorers'], 'dir')

      mkdir([wDir, filesep, 'Scorers']);   
   end
   
   nPhantoms = size(phantoms,2);
   nIons = size(Ions,2);
   
   for k=1:nPhantoms
      currDir = [wDir, filesep, 'Scorers', filesep, phantoms(k).name];
      name = phantoms(k).name;
      XBins = ceil(phantoms(k).dimension(2)/phantoms(k).resolution(2));
      YBins = ceil(phantoms(k).dimension(1)/phantoms(k).resolution(1));
      ZBins = ceil(phantoms(k).dimension(3)/phantoms(k).resolution(3));
      
      if ~exist(currDir, 'dir')
         mkdir(currDir);   
      end
      
      %cd(currDir);
      %Write geometryPhantom
      fileName = [currDir,filesep,'Geometry_scorer_', name, '.txt'];
      fID = fopen(fileName, 'w');
      fprintf(fID, 's:Ge/%s/Parent="World" \n',name);
      fprintf(fID, 's:Ge/%s/Type = "TsBox" \n',name);
      fprintf(fID, 's:Ge/%s/Material = "G4_WATER" \n',name);
      fprintf(fID, 'd:Ge/%s/HLX = %f mm \n',name,phantoms(k).dimension(2)/2);
      fprintf(fID, 'd:Ge/%s/HLY = %f mm \n',name,phantoms(k).dimension(1)/2);
      fprintf(fID, 'd:Ge/%s/HLZ = %f mm \n',name,phantoms(k).dimension(3)/2);
      fprintf(fID, 'i:Ge/%s/XBins = %d \n',name,XBins);
      fprintf(fID, 'i:Ge/%s/YBins = %d \n',name,YBins);
      fprintf(fID, 'i:Ge/%s/ZBins = %d \n',name,ZBins);
      fprintf(fID, 'd:Ge/%s/TransX = %f mm \n',name,phantoms(k).resolution(2)/2);
      fprintf(fID, 'd:Ge/%s/TransY = %f mm \n',name,(phantoms(k).dimension(1)/2 + phantoms(1).positionDepth(k)));
      fprintf(fID, 'd:Ge/%s/TransZ = %f mm \n',name,phantoms(k).resolution(3)/2);
      fprintf(fID, 'b:Ge/%s/IsParallel = "True" \n',name);
      fprintf(fID, 'includeFile = %s', ['./beamSetup_', topas_cfg.label, '_field1.txt']);
      fprintf(fID, '\n');
      fclose(fID);
      
      fileName = [currDir,filesep,'scorer_', name, '.txt'];
      fID = fopen(fileName, 'w');
      if includePDD
        fprintf(fID, 'includeFile = %s \n', ['./Scorers/', phantoms(k).name,'/PDD_', name, '.txt']);
      end

      if includeDS
          for m=1:nIons
             fprintf(fID, 'includeFile = %s \n', ['./Scorers/',phantoms(k).name, '/S_', name, '_', Ions{m}, '.txt']);
          end
      end
      
      if includeSP
          for m=1:nIons
            fprintf(fID, 'includeFile = %s \n', ['./Scorers/',phantoms(k).name, '/SP_', name, '_', Ions{m}, '.txt']);
          end
      end

      if includeEdEventSpectrum
          for m=1:nIons
            fprintf(fID, 'includeFile = %s \n', ['./Scorers/',phantoms(k).name, '/SPeD_', name, '_', Ions{m}, '.txt']);

          end
      end

      if includeProtonLET && strcmp([Ions{:}],'protons')
        fprintf(fID, 'includeFile = %s \n', ['./Scorers/', phantoms(k).name,'/ProtonLET_', name, '.txt']);
      end
      %fprintf(fID, '#includeFile = %s', ['./beamSetup_', topas_cfg.label, '_field1.txt']);
      fclose(fID);

      if includePDD
        if ~exist(strcat(topas_cfg.workingDir,filesep, 'Output\PDD',filesep, name), 'dir')
          mkdir(strcat(topas_cfg.workingDir,filesep, 'Output\PDD',filesep, name));
        end
        writePDDscorer(topas_cfg, phantoms(k),currDir);
      end

      if includeDS
          for m=1:nIons
            if ~exist(strcat(topas_cfg.workingDir,filesep, 'Output\DS',filesep, name,filesep, strcat('Ion_',num2str(m))), 'dir')
                mkdir(strcat(topas_cfg.workingDir,filesep, 'Output\DS',filesep, name,filesep, strcat('Ion_',num2str(m))));
            end
          end
        writeDSscorers(EParam, Ions, topas_cfg, phantoms(k),currDir);
      end

       if includeEdEventSpectrum
          for m=1:nIons
            if ~exist(strcat(topas_cfg.workingDir,filesep, 'Output\SPeD',filesep, name,filesep, strcat('Ion_',num2str(m))), 'dir')
                mkdir(strcat(topas_cfg.workingDir,filesep, 'Output\SPeD',filesep, name,filesep, strcat('Ion_',num2str(m))));
            end
          end
        writeSPeDscorers(EParam, Ions, topas_cfg, phantoms(k),currDir);
      end

      if includeSP
         for m=1:nIons
            if ~exist(strcat(topas_cfg.workingDir,filesep, 'Output\SP',filesep, name,filesep, strcat('Ion_',num2str(m))), 'dir')
                mkdir(strcat(topas_cfg.workingDir,filesep, 'Output\SP',filesep, name,filesep, strcat('Ion_',num2str(m))));
            end
          end
        writeSPScorers(EParam, Ions, topas_cfg, phantoms(k),currDir);
      end

      if includeProtonLET && (strcmp([Ions{:}], 'protons'))
        if ~exist(strcat(topas_cfg.workingDir,filesep, 'Output\ProtonLET',filesep, name), 'dir')
          mkdir(strcat(topas_cfg.workingDir,filesep, 'Output\ProtonLET',filesep, name));
        end
        writeProtonLETScorer(topas_cfg, phantoms(k),currDir);
      end
  end
   
end