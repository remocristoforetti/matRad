function WriteProtonLUT(LUTDir,wDir)
    load(LUTDir);
    Idx = find(E>0.1);
    E = E(Idx);
    G = G(Idx);
    
    fid = fopen(strcat(wDir, filesep, 'LUTProtons.txt'), 'w');
    fprintf(fid, 'dv:Sc/Patient/ZMIX/KineticEnergyPerNucleon = %i ', size(E,2));
    fprintf(fid, num2str(E));
    fprintf(fid,' MeV');
    fprintf(fid, '\n');
    fprintf(fid, 'dv:Sc/Patient/ZMIX/Zd = %i ', size(G,2));
    fprintf(fid, num2str(G));
    fprintf(fid, ' Gy');
    fprintf(fid, '\n');
    
    fclose(fid);

end