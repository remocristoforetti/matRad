function writeScoringTreeDirectory(phantoms, Ions, topas_cfg, EParam)

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
      fprintf(fID, '#s:Ge/%s/Material = G4_WATER \n',name);
      fprintf(fID, 'd:Ge/%s/HLX = %f mm \n',name,phantoms(k).dimension(2)/2);
      fprintf(fID, 'd:Ge/%s/HLY = %f mm \n',name,phantoms(k).dimension(1)/2);
      fprintf(fID, 'd:Ge/%s/HLZ = %f mm \n',name,phantoms(k).dimension(3)/2);
      fprintf(fID, 'i:Ge/%s/XBins = %d \n',name,XBins);
      fprintf(fID, 'i:Ge/%s/YBins = %d \n',name,YBins);
      fprintf(fID, 'i:Ge/%s/ZBins = %d \n',name,ZBins);
      fprintf(fID, 'd:Ge/%s/TransX = %f mm \n',name,phantoms(k).resolution(2)/2);
      fprintf(fID, 'd:Ge/%s/TransY = %f mm \n',name,-(phantoms(k).dimension(1)/2 + phantoms(1).positionDepth(k)));
      fprintf(fID, 'd:Ge/%s/TransZ = %f mm \n',name,phantoms(k).resolution(3)/2);
      fprintf(fID, 'b:Ge/%s/IsParallel = "True" \n',name);
      fprintf(fID, 'includeFile = %s', ['./beamSetup_', topas_cfg.label, '_field1.txt']);
      fprintf(fID, '\n');
      fclose(fID);
      
      fileName = [currDir,filesep,'scorer_', name, '.txt'];
      fopen(fileName, 'w');
      fprintf(fID, 'includeFile = %s \n', ['./Scorers/', phantoms(k).name,'/PDD_', name, '.txt']);
      for m=1:nIons
         fprintf(fID, 'includeFile = %s \n', ['./Scorers/',phantoms(k).name, '/S_', name, '_', Ions{m}, '.txt']);
      end
      fprintf(fID, '#includeFile = %s', ['./beamSetup_', topas_cfg.label, '_field1.txt']);
      fclose(fID);
      writePDDscorer(topas_cfg, phantoms(k),currDir);
      writeDSscorers(EParam, Ions, topas_cfg, phantoms(k),currDir);
   end
   
end