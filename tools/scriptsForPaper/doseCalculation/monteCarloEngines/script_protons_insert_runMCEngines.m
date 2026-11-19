%% Script: Run the three MC codes (TOPAS, MCsquare, FRED) on pre-written input files
%% set matRad runtime configuration
clear;
matRad_rc;
matRad_cfg = MatRad_Config.instance();

mainDir = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','doseCalculation','MCengines','protonsInsert');

% choose which codes to run
runTopas    = true;
runMCsquare = true;
runFred     = true;
fredUseGPU  = true; % set true if FRED has a GPU available

startDir = pwd;
% cleanupObj = onCleanup(@() cd(startDir)); % always return to the starting folder

%% TOPAS
if runTopas
    try
        topasDir = fullfile(mainDir,'TOPAS');
        topasPlanFiles = dir(fullfile(topasDir,'matRad_plan_field*_run*.txt'));
        if isempty(topasPlanFiles)
            matRad_cfg.dispError('No TOPAS plan files (matRad_plan_field*_run*.txt) found in %s',topasDir);
        end

        engineTopas = DoseEngines.matRad_TopasMCEngine(); % only used for its default exec command
        engineTopas.topasExecCommand = 'wsl export TOPAS_G4_DATA_DIR=~/Applications/GEANT4/G4DATA; ~/Applications/TOPAS/OpenTOPAS-install/bin/topas';
        cd(topasDir);
        for i = 1:numel(topasPlanFiles)
            [~,fname] = fileparts(topasPlanFiles(i).name);
            topasCall = sprintf('%s %s.txt > %s.out > %s.log',engineTopas.topasExecCommand,fname,fname,fname);
            status = system(topasCall,'-echo');
            if status ~= 0
                matRad_cfg.dispError('TOPAS exited with error code %d for %s',status,fname);
            end
        end
        cd(startDir);
    catch ME
        cd(startDir);
        matRad_cfg.dispWarning('TOPAS run failed: %s',ME.message);
    end
end

%% MCsquare
if runMCsquare
    try
        mcSquareDir = fullfile(mainDir,'MCsquare');
        mcSquareConfigFile = fullfile(mcSquareDir,'MCsquareConfig.txt');
        if ~isfile(mcSquareConfigFile)
            matRad_cfg.dispError('No MCsquare config file found in %s',mcSquareDir);
        end

        mcSquareBinFolder = fullfile(matRad_cfg.matRadRoot,'thirdParty','MCsquare','bin');
        cd(mcSquareBinFolder);
        if ~isfile(fullfile(mcSquareBinFolder,'Materials','list.dat'))
            unzip('Materials.zip');
        end
        [binaryFound,mcSquareBinary] = DoseEngines.matRad_ParticleMCsquareEngine.checkBinaries();

        % MCsquare stores file paths in fixed-length buffers and crashes
        % silently (exit code 0xC0000409, no output) when a path in the
        % config is longer than ~125 characters. matRad writes absolute
        % paths, so write a copy of the config with all paths made
        % relative to the bin folder (the call is run from there).
        rootFwd    = strrep(matRad_cfg.matRadRoot,'\','/');
        relRoot    = '../../..'; % bin -> thirdParty/MCsquare/bin -> matRadRoot
        mcSquareDirFwd = strrep(mcSquareDir,'\','/');
        if ~startsWith(mcSquareDirFwd,rootFwd)
            matRad_cfg.dispError('MCsquare folder must be inside the matRad root (%s) to use relative paths',rootFwd);
        end
        configText = fileread(mcSquareConfigFile);
        configText = strrep(configText,rootFwd,relRoot);
        runConfigFile = fullfile(mcSquareDir,'MCsquareConfig_run.txt');
        fid = fopen(runConfigFile,'w');
        fwrite(fid,configText);
        fclose(fid);
        runConfigRel = [relRoot mcSquareDirFwd(numel(rootFwd)+1:end) '/MCsquareConfig_run.txt'];

        mcSquareCall = sprintf('%s "%s"',mcSquareBinary,runConfigRel);
        status = system(mcSquareCall,'-echo');
        cd(startDir);
        if status ~= 0
            matRad_cfg.dispError('MCsquare exited with error code %d',status);
        end
        if ~isfile(fullfile(mcSquareDir,'output','Dose.mhd'))
            matRad_cfg.dispError('MCsquare finished but wrote no Dose.mhd to %s',fullfile(mcSquareDir,'output'));
        end
    catch ME
        cd(startDir);
        matRad_cfg.dispWarning('MCsquare run failed: %s',ME.message);
    end
end

%% FRED
if runFred
    try
        fredRunDir = fullfile(mainDir,'FRED','MCrun');
        if ~isfile(fullfile(fredRunDir,'fred.inp'))
            matRad_cfg.dispError('No fred.inp found in %s',fredRunDir);
        end

        flags = '-f fred.inp';
        if ~fredUseGPU
            flags = ['-nogpu ' flags];
        end
        fredCall = [strtrim(DoseEngines.matRad_ParticleFREDEngine.cmdCall) ' ' flags];

        cd(fredRunDir);
        status = system(fredCall,'-echo');
        cd(startDir);
        if status ~= 0
            matRad_cfg.dispError('FRED exited with error code %d',status);
        end
        matRad_cfg.dispInfo('FRED finished, results in %s\n',fullfile(fredRunDir,'out'));
    catch ME
        cd(startDir);
        matRad_cfg.dispWarning('FRED run failed: %s',ME.message);
    end
end
