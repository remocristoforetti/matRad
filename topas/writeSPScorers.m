function writeSPScorers(EParam, Ions, topas_cfg, phantom, currDir)
   nIons = size(Ions,2);
   IonsZ = [];
   for k=1:nIons
      
      switch Ions{k}
         case 'protons'
            IonsZ =  1;
            A = 1;
         case 'He'
            IonsZ =  2;
            A = 4;
         case 'Li'
            IonsZ =  3;
            A = 7;
         case 'Be'
            IonsZ =  4;
            A = 9;
         case 'B'
            IonsZ =  5;
            A = 11;
         case 'C'
            IonsZ =  6;
            A = 12;
      end

      fileName = [currDir, filesep,'SP_',phantom.name, '_', Ions{k}, '.txt'];
      fID = fopen(fileName, 'w');
      
      fprintf(fID,'#-- E Spectrum Scorer for %s, this scoreer scores the particle fluence, not the event fluence \n',Ions{k});
      fprintf(fID,'s:Sc/%s/%s/%s/Quantity = "Fluence" \n',phantom.name,'SP',Ions{k});
      fprintf(fID,'s:Sc/%s/%s/%s/OutputType = "csv" \n',phantom.name,'SP',Ions{k});
      fprintf(fID,'s:Sc/%s/%s/%s/Component = "%s" \n',phantom.name,'SP',Ions{k},phantom.name);
      fprintf(fID,'s:Sc/%s/%s/%s/OutputFile = "./Output/SP/%s/Ion_%1d/SP_%s" \n',phantom.name,'SP',Ions{k},phantom.name,IonsZ,Ions{k});
      fprintf(fID,'i:Sc/%s/%s/%s/OnlyIncludeParticlesOfAtomicNumber = %1d \n',phantom.name,'SP',Ions{k},IonsZ);
      fprintf(fID,'s:Sc/%s/%s/%s/IfOutputFileAlreadyExists = "Overwrite" \n',phantom.name,'SP',Ions{k});
      fprintf(fID,'b:Sc/%s/%s/%s/OutputAfterRun = "True" \n',phantom.name,'SP',Ions{k});
      fprintf(fID, '\n');
      
      fprintf(fID,'d:Sc/%s/%s/%s/EBinMax = %3.4f MeV \n',phantom.name,'SP',Ions{k},EParam(k).EMax*A);
      fprintf(fID,'d:Sc/%s/%s/%s/EBinMin = %3.4f MeV \n',phantom.name,'SP',Ions{k},EParam(k).EMin*A);
      fprintf(fID,'i:Sc/%s/%s/%s/EBins = %4u \n',phantom.name,'SP',Ions{k},EParam(k).nEBins);
      fprintf(fID, '\n');
      
      fprintf(fID,'i:Sc/%s/%s/%s/XBins = 1\n',phantom.name,'SP',Ions{k});
      fprintf(fID,'i:Sc/%s/%s/%s/YBins = Ge/%s/YBins \n',phantom.name,'SP',Ions{k},phantom.name);
      fprintf(fID,'i:Sc/%s/%s/%s/ZBins = 1 \n',phantom.name,'SP',Ions{k});
      fprintf(fID,'includeFile = %s',['./Scorers/', phantom.name, '/Geometry_scorer_', phantom.name, '.txt']);
      fclose(fID);
   end
end