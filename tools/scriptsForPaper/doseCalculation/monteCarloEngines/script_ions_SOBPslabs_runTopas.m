%% Script: Run TOPAS on the helium/carbon/oxygen SOBP input files (entrance slabs)
%% set matRad runtime configuration
clear;
matRad_rc;
matRad_cfg = MatRad_Config.instance();

ions = {'helium','carbon','oxygen'};

topasExecCommand = 'wsl export TOPAS_G4_DATA_DIR=~/Applications/GEANT4/G4DATA; ~/Applications/TOPAS/OpenTOPAS-install/bin/topas';

startDir = pwd;

for i = 1:numel(ions)
    try
        topasDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','doseCalculation','MCengines',[ions{i} 'SOBPslabs'],'TOPAS');
        topasPlanFiles = dir(fullfile(topasDir,'matRad_plan_field*_run*.txt'));
        if isempty(topasPlanFiles)
            matRad_cfg.dispError('No TOPAS plan files (matRad_plan_field*_run*.txt) found in %s',topasDir);
        end

        cd(topasDir);
        for j = 1:numel(topasPlanFiles)
            [~,fname] = fileparts(topasPlanFiles(j).name);
            topasCall = sprintf('%s %s.txt > %s.out > %s.log',topasExecCommand,fname,fname,fname);
            matRad_cfg.dispInfo('%s: calling TOPAS: %s\n',ions{i},topasCall);
            status = system(topasCall,'-echo');
            if status ~= 0
                matRad_cfg.dispError('TOPAS exited with error code %d for %s',status,fname);
            end
        end
        cd(startDir);
    catch ME
        cd(startDir);
        matRad_cfg.dispWarning('%s: TOPAS run failed: %s',ions{i},ME.message);
    end
end
