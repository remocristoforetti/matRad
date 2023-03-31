function writeScoringTreeDirectory_singlePhantom(phantom, ions, topas_cfg)
   wDir = topas_cfg.workingDir;
   %make Scorer directory
   if ~exist([wDir,filesep,'Scorers'], 'dir')
      mkdir([wDir, filesep, 'Scorers']);   

   end
   
   nIons = size(ions,2);
   currDir = [wDir, filesep, 'Scorers',filesep, phantom{1,5}.name];
   name = phantom{1,5}.name;
   
   XBins = ceil(phantom{1,5}.dimension.x/phantom{1,5}.resolution.x);
   YBins = ceil(phantom{1,5}.dimension.y/phantom{1,5}.resolution.y);
   ZBins = ceil(phantom{1,5}.dimension.z/phantom{1,5}.resolution.z);
      
  %Write geometryPhantom
  fileName = [currDir,filesep,'Geometry_scorer_', name, '.txt'];

  if ~exist(currDir, 'dir')
      mkdir(currDir);   
   end
  fID = fopen(fileName, 'w');
  fprintf(fID, 's:Ge/%s/Parent="World" \n',name);
  fprintf(fID, 's:Ge/%s/Type = "TsBox" \n',name);
  fprintf(fID, 's:Ge/%s/Material = "G4_WATER" \n',name);
  fprintf(fID, 'd:Ge/%s/HLX = %f mm \n',name,phantom{1,5}.dimension.x/2);
  fprintf(fID, 'd:Ge/%s/HLY = %f mm \n',name,phantom{1,5}.dimension.y/2);
  fprintf(fID, 'd:Ge/%s/HLZ = %f mm \n',name,phantom{1,5}.dimension.z/2);
  fprintf(fID, 'i:Ge/%s/XBins = %d \n',name,XBins);
  fprintf(fID, 'i:Ge/%s/YBins = %d \n',name,YBins);
  fprintf(fID, 'i:Ge/%s/ZBins = %d \n',name,ZBins);
  fprintf(fID, 'd:Ge/%s/TransX = %f mm \n',name,phantom{1,5}.resolution.x/2);
  fprintf(fID, 'd:Ge/%s/TransY = %f mm \n',name,phantom{1,5}.dimension.y/2);
  fprintf(fID, 'd:Ge/%s/TransZ = %f mm \n',name,phantom{1,5}.resolution.z/2);
  fprintf(fID, 'b:Ge/%s/IsParallel = "True" \n',name);
  fprintf(fID, 'includeFile = %s', ['./beamSetup_', topas_cfg.label, '_field1.txt']);
  fprintf(fID, '\n');
  fclose(fID);
      
   fileName = [currDir,filesep,'scorer_', name, '.txt'];
   fID = fopen(fileName, 'w');
   if phantom{1,5}.includePDD
        fprintf(fID, 'includeFile = %s \n', ['./Scorers/', name, '/PDD_', name, '.txt']);
   end

  if phantom{1,5}.includeDS
      for m=1:nIons
         fprintf(fID, 'includeFile = %s \n', ['./Scorers/', name,'/S_', name, '_', ions{m}, '.txt']);
      end
  end
      
  if phantom{1,5}.includeSP
      for m=1:nIons
        fprintf(fID, 'includeFile = %s \n', ['./Scorers/', name,'/SP_', name, '_', ions{m}, '.txt']);

      end
  end

  if phantom{1,5}.includeEdEventSpectrum
      for m=1:nIons
        fprintf(fID, 'includeFile = %s \n', ['./Scorers/', name,'/SPeD_', name, '_', ions{m}, '.txt']);

      end
  end

  if phantom{1,5}.includeProtonLET && strcmp([ions{:}],'protons')
    fprintf(fID, 'includeFile = %s \n', ['./Scorers/', name,'/ProtonLET_', name, '.txt']);
  end

 if phantom{1,5}.includeSPdose
      for m=1:nIons
        fprintf(fID, 'includeFile = %s \n', ['./Scorers/', name,'/SPdose_', name, '_', ions{m}, '.txt']);

      end
 
 end

 if phantom{1,5}.includeEventsCounter
      for m=1:nIons
        fprintf(fID, 'includeFile = %s \n', ['./Scorers/', name,'/Phi_', name, '_', ions{m}, '.txt']);

      end
  end

  %fprintf(fID, '#includeFile = %s', ['./beamSetup_', topas_cfg.label, '_field1.txt']);
  fclose(fID);

  if phantom{1,5}.includePDD
    if ~exist(strcat(topas_cfg.workingDir,filesep, 'Output\PDD',filesep, name), 'dir')
      mkdir(strcat(topas_cfg.workingDir,filesep, 'Output\PDD',filesep, name));
    end
    writePDDscorer(topas_cfg, phantom{1,5},currDir);
  end

  if phantom{1,5}.includeDS
      for m=1:nIons
        if ~exist(strcat(topas_cfg.workingDir,filesep, 'Output\DS',filesep, name,filesep, strcat('Ion_',num2str(m))), 'dir')
            mkdir(strcat(topas_cfg.workingDir,filesep, 'Output\DS',filesep, name,filesep, strcat('Ion_',num2str(m))));
        end
      end
    writeDSscorers(phantom{1,4}, ions, topas_cfg, phantom{1,5},currDir);
  end

  if phantom{1,5}.includeEdEventSpectrum
      for m=1:nIons
        if ~exist(strcat(topas_cfg.workingDir,filesep, 'Output\SPeD',filesep, name,filesep, strcat('Ion_',num2str(m))), 'dir')
            mkdir(strcat(topas_cfg.workingDir,filesep, 'Output\SPeD',filesep, name,filesep, strcat('Ion_',num2str(m))));
        end
      end
    writeSPeDscorers(phantom{1,4}, ions, topas_cfg, phantom{1,5},currDir);
  end

 
 if phantom{1,5}.includeSP
     for m=1:nIons
        if ~exist(strcat(topas_cfg.workingDir,filesep, 'Output\SP',filesep, name,filesep, strcat('Ion_',num2str(m))), 'dir')
            mkdir(strcat(topas_cfg.workingDir,filesep, 'Output\SP',filesep, name,filesep, strcat('Ion_',num2str(m))));
        end
      end
    writeSPScorers(phantom{1,4}, ions, topas_cfg, phantom{1,5},currDir);
 end

  if phantom{1,5}.includeProtonLET && (strcmp([ions{:}], 'protons'))
    if ~exist(strcat(topas_cfg.workingDir,filesep, 'Output\ProtonLET',filesep, name), 'dir')
      mkdir(strcat(topas_cfg.workingDir,filesep, 'Output\ProtonLET',filesep, name));
    end
    writeProtonLETScorer(topas_cfg, phantom{1,5},currDir);
  end

  if phantom{1,5}.includeSPdose
     for m=1:nIons
        if ~exist(strcat(topas_cfg.workingDir,filesep, 'Output\SPdose',filesep, name,filesep, strcat('Ion_',num2str(m))), 'dir')
            mkdir(strcat(topas_cfg.workingDir,filesep, 'Output\SPdose',filesep, name,filesep, strcat('Ion_',num2str(m))));
        end
      end
    writeSPdoseScorers(phantom{1,4}, ions, topas_cfg, phantom{1,5},currDir);
  end

  if phantom{1,5}.includeEventsCounter
      for m=1:nIons
        if ~exist(strcat(topas_cfg.workingDir,filesep, 'Output\Phi',filesep, name,filesep, strcat('Ion_',num2str(m))), 'dir')
            mkdir(strcat(topas_cfg.workingDir,filesep, 'Output\Phi',filesep, name,filesep, strcat('Ion_',num2str(m))));
        end
      end
    writeEventsCounterScorers(phantom{1,4}, ions, topas_cfg, phantom{1,5},currDir);
  end
end