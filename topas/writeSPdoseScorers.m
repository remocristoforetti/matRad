function writeSPdoseScorers(EParam, Ions, topas_cfg, phantom, currDir)

   nIons = size(Ions,2);
   IonsZ = [];
   for k=1:nIons
      
      switch Ions{k}
         case 'protons'
            IonsZ =  [IonsZ, 1];
         case 'He'
            IonsZ =  [IonsZ, 2];   
         case 'Li'
            IonsZ =  [IonsZ, 3];
         case 'Be'
            IonsZ =  [IonsZ, 4];
         case 'B'
            IonsZ =  [IonsZ, 5];
         case 'C'
            IonsZ =  [IonsZ, 6];
      end
   end
   
   for k=1:nIons
      fileName = [currDir, filesep,'SPdose_',phantom.name, '_', Ions{k}, '.txt'];
      fID = fopen(fileName, 'w');
      
      fprintf(fID,'#-- Ed Scorer for %s, this scores the energy deposioted by all the events within the voxel \n',Ions{k});
      %fprintf(fID,'s:Sc/%s/%s/%s/Quantity = "TsEdEventSpectrum" \n',phantom.name,'SPeD',Ions{k});
      fprintf(fID,'s:Sc/%s/%s/%s/Quantity = "TsDoseEventsCounter" \n',phantom.name,'SPdose',Ions{k});
      fprintf(fID,'s:Sc/%s/%s/%s/OutputType = "DICOM" \n',phantom.name,'SPdose',Ions{k});
      fprintf(fID,'s:Sc/%s/%s/%s/Component = "%s" \n',phantom.name,'SPdose',Ions{k},phantom.name);
      fprintf(fID,'s:Sc/%s/%s/%s/OutputFile = "./Output/SPdose/%s/Ion_%1d/SPdose_%s" \n',phantom.name,'SPdose',Ions{k},phantom.name,IonsZ(k),Ions{k});
      fprintf(fID,'i:Sc/%s/%s/%s/OnlyIncludeParticlesOfAtomicNumber = %1d \n',phantom.name,'SPdose',Ions{k},IonsZ(k));
      fprintf(fID,'s:Sc/%s/%s/%s/IfOutputFileAlreadyExists = "Overwrite" \n',phantom.name,'SPdose',Ions{k});
      fprintf(fID,'b:Sc/%s/%s/%s/OutputAfterRun = "True" \n',phantom.name,'SPdose',Ions{k});
      
      %"ReferencedSubScorer_" is fixed, "nEvents_subScorer" is reference to
      %subscorer in extensions (third parameter of InstantiateSubScorer)
      %and then need full name of scorer -> In this case is the events
      %counter (instantiated in other parameter file)
      %fprintf(fID,'s:Sc/%s/%s/%s/ReferencedSubScorer_nEvents_subScorer = "%s/%s/%s" \n',phantom.name,'SPdose',Ions{k},phantom.name,'Phi',Ions{k});

      fprintf(fID, '\n');     
      
      fprintf(fID,'d:Sc/%s/%s/%s/Emax = %3.4f MeV \n',phantom.name,'SPdose',Ions{k},EParam(k).EMax);
      fprintf(fID,'d:Sc/%s/%s/%s/Emin = %3.4f MeV \n',phantom.name,'SPdose',Ions{k},EParam(k).EMin);
      fprintf(fID,'i:Sc/%s/%s/%s/nEBins = %4u \n',phantom.name,'SPdose',Ions{k},EParam(k).nEBins);
      fprintf(fID, '\n');
      
%       fprintf(fID,'d:Sc/%s/%s/%s/Edepmax = %3.4f MeV \n',phantom.name, 'SPeD', Ions{k},EParam(k).EdMax);
%       fprintf(fID,'d:Sc/%s/%s/%s/Edepmin = %3.4f MeV \n',phantom.name, 'SPeD',Ions{k},EParam(k).EdMin);
%       fprintf(fID,'i:Sc/%s/%s/%s/nEdepBins = %4u \n',phantom.name, 'SPeD',Ions{k},EParam(k).nEdBins);
%       fprintf(fID, '\n');
      
      fprintf(fID,'i:Sc/%s/%s/%s/XBins = 1 \n',phantom.name,'SPdose',Ions{k});
      fprintf(fID,'i:Sc/%s/%s/%s/YBins = Ge/%s/YBins \n',phantom.name,'SPdose',Ions{k},phantom.name);
      fprintf(fID,'i:Sc/%s/%s/%s/ZBins = Sc/%s/%s/%s/nEBins + 1 \n',phantom.name,'SPdose',Ions{k},phantom.name,'SPdose',Ions{k});
      fprintf(fID, 'includeFile = %s',['./Scorers/', phantom.name, '/Geometry_scorer_', phantom.name, '.txt']);
      fclose(fID);
   end
end