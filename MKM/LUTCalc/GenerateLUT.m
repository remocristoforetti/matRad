function time = GenerateLUT(domain, Integral, track, Z, E, wDir)
    for ionZ = Z
        tic
        track.Z_ion = ionZ;
        zD = [];
        reverseStr = '';

        for k=1:length(E)

            track.Ek = E(k)*track.A;
            Integral.computeZD(track,domain);
            zD = [zD, Integral.gamma];
            
            percentDone = 100 * k / length(E);
            msg = sprintf('Computing LUT for Z = %d: %3.1f', ionZ, percentDone);
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));        
        
        end
        time = toc;
        save(strcat(wDir,filesep, 'LUT_Z_', num2str(ionZ), '.mat'), 'E', 'zD');
        fprintf('\n');

   end
end