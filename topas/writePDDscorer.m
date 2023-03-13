function writePDDscorer(topas_cfg, phantom,currDir)
   fileName = [currDir, filesep, 'PDD_', phantom.name, '.txt'];
   fID = fopen(fileName, 'w');
   fprintf(fID, 's:Sc/%s/DoseToMedium/Quantity = "DoseToMedium" \n', phantom.name);
   fprintf(fID, 's:Sc/%s/DoseToMedium/OutputType = "DICOM" \n', phantom.name);

   fprintf(fID, 's:Sc/%s/DoseToMedium/Component = "%s" \n', phantom.name, phantom.name);
   fprintf(fID, 'sv:Sc/%s/DoseToMedium/Report = 1 "Sum" \n', phantom.name);
   fprintf(fID, 's:Sc/%s/DoseToMedium/IfOutputFileAlreadyExists = "Overwrite" \n', phantom.name);
   fprintf(fID, 'b:Sc/%s/DoseToMedium/OutputToConsole = "False" \n', phantom.name);
   fprintf(fID, 's:Sc/%s/DoseToMedium/OutputFile = "./Output/PDD/%s/physicalDose" \n', phantom.name, phantom.name);
   fprintf(fID, 'b:Sc/%s/DoseToMedium/OutputAfterRun = "True" \n', phantom.name);
   fprintf(fID, 'includeFile = ./Scorers/%s/Geometry_scorer_%s.txt',phantom.name,phantom.name);
   fclose(fID);
end