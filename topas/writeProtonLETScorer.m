function writeProtonLETScorer(topas_cfg, phantom,currDir)
   fileName = [currDir, filesep, 'ProtonLET_', phantom.name, '.txt'];
   fID = fopen(fileName, 'w');
   fprintf(fID, 's:Sc/%s/ProtonLET/Quantity = "ProtonLET" \n', phantom.name);
   fprintf(fID, 's:Sc/%s/ProtonLET/OutputType = "DICOM" \n', phantom.name);

   fprintf(fID, 's:Sc/%s/ProtonLET/Component = "%s" \n', phantom.name, phantom.name);
   fprintf(fID, 'sv:Sc/%s/ProtonLET/Report = 1 "Sum" \n', phantom.name);
   fprintf(fID, 's:Sc/%s/ProtonLET/IfOutputFileAlreadyExists = "Overwrite" \n', phantom.name);
   fprintf(fID, 'b:Sc/%s/ProtonLET/OutputToConsole = "False" \n', phantom.name);
   fprintf(fID, 's:Sc/%s/ProtonLET/OutputFile = "./Output/ProtonLET/%s/protonLET" \n', phantom.name, phantom.name);
   fprintf(fID, 'b:Sc/%s/ProtonLET/OutputAfterRun = "True" \n', phantom.name);

   fprintf(fID,'i:Sc/%s/%s/XBins = 1 \n',phantom.name,'ProtonLET');
   fprintf(fID,'i:Sc/%s/%s/YBins = Ge/%s/YBins \n',phantom.name,'ProtonLET',phantom.name);
   
   fprintf(fID,'i:Sc/%s/%s/ZBins = 1 \n',phantom.name,'ProtonLET');
   
   fprintf(fID, 'includeFile = ./Scorers/%s/Geometry_scorer_%s.txt',phantom.name,phantom.name);
   fclose(fID);
end