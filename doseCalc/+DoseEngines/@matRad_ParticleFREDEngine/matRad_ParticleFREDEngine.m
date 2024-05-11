classdef matRad_ParticleFREDEngine < DoseEngines.matRad_MonteCarloEngineAbstract
% Engine for particle dose calculation using monte carlo calculation
% specificly the FRED MC code ()
% for more informations see superclass
% DoseEngines.matRad_MonteCarloEngineAbstract
%
% pln.propDoseCalc fields:
% 
% HUclamping:           [b] allows for clamping of HU table
% HUtable:              [s] 'internal', 'matRad_water', 'matRad_water_78'
% exportCalculation     [b] t:Only write simulation paramter files, f: run
%                             FRED
% sourceModel           [s] see AvailableSourceModels, {'gaussian', 'emittance', 'sigmaSqrModel'}
% useWSL                [b]
% useGPU                [b]
% useWaterPhantom       [b] overwrite ct phantom with uniform water phantom
%                           (same resolution as doseGrid)
% roomMaterial          [s] vacuum, Air


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    properties (Constant)
        
        possibleRadiationModes = {'protons', 'carbon'};
        name = 'FRED';
        shortName = 'FRED';
    end

    properties (SetAccess = protected, GetAccess = public)
        
        currFolder = pwd; %folder path when set
        
        %nbThreads; %number of threads for MCsquare, 0 is all available

        %useInternalHUConversion;
        noozleToAxis;
        scorers = {'Dose'};
     
       
        HUtable;
        defaultHUtable = 'matRad_water.txt';

        availableRBEmodels = {'MCN'};
        calcBioDose;
        vAlphaX;
        vBetaX;
        vABratio;

        sourceModel;
        AvailableSourceModels = {'gaussian', 'emittance', 'sigmaSqrModel'};

        useWaterPhantom;
        roomMaterial;
        useWSL = false;
        useGPU;
        currentVersion;
        availableVersions;
        fragDataLib;
    end

    properties
        exportCalculation = false;
        calcLET;
        %RBEmodel = 'none';
        constantRBE = NaN;              % constant RBE value
        HUclamping = false;
    end

    properties (SetAccess = private, Hidden)
        patientFilename      = 'CTpatient.mhd';
        runInputFilename     = 'fred.inp';
        
        regionsFilename      = 'regions.inp';
        funcsFilename        = 'funcs.inp';
        planFilename         = 'plan.inp';
        fieldsFilename       = 'fields.inp';
        layersFilename       = 'layers.inp';
        beamletsFilename     = 'beamlets.inp';
        planDeliveryFilename = 'planDelivery.inp';

        hLutLimits = [-1000,1375];
        
        conversionFactor = 1e6;
        %conversionFactor = 1;
        planDeliveryTemplate = 'planDelivery.txt';

        FREDrootFolder;

        MCrunFolder;
        inputFolder;
        regionsFolder;
        planFolder;

        cmdCall;

    end
    
    methods
        
        function this = matRad_ParticleFREDEngine(pln)
            % Constructor
            %
            % call
            %   engine = DoseEngines.matRad_DoseEngineFRED(ct,stf,pln,cst)
            %
            % input
            %   ct:                         matRad ct struct
            %   stf:                        matRad steering information struct
            %   pln:                        matRad plan meta information struct
            %   cst:                        matRad cst struct
            
            matRad_cfg = MatRad_Config.instance();
            if nargin < 1
                pln = [];
            end

            % call superclass constructor
            this = this@DoseEngines.matRad_MonteCarloEngineAbstract(pln);


            if nargin > 0
                if isfield(pln,'propDoseCalc') && isfield(pln.propDoseCalc,'HUclamping')
                    this.HUclamping = pln.propDoseCalc.HUclamping;
                end

                if isfield(pln,'propDoseCalc') && isfield(pln.propDoseCalc,'HUtable')
                    this.HUtable = pln.propDoseCalc.HUtable;
                else
                    this.HUtable = this.defaultHUtable;
                end

                if isfield(pln,'propDoseCalc') && isfield(pln.propDoseCalc,'exportCalculation')
                    this.exportCalculation = pln.propDoseCalc.exportCalculation; 
                else
                    this.exportCalculation = matRad_cfg.propMC.defaultExternalCalculation;
                end

                if isfield(pln,'propDoseCalc') && isfield(pln.propDoseCalc,'sourceModel')
                    this.sourceModel = pln.propDoseCalc.sourceModel; 
                else
                    this.sourceModel = matRad_cfg.propMC.defaultSourceModel;
                end

                if isfield(pln, 'propDoseCalc') && isfield(pln.propDoseCalc,'useWSL')
                    this.useWSL = pln.propDoseCalc.useWSL;
                else
                    this.useWSL = false;
                end

                if isfield(pln, 'propDoseCalc') && isfield(pln.propDoseCalc,'useGPU')
                    this.useGPU = pln.propDoseCalc.useGPU;
                else
                    this.useGPU = true;
                end

                if isfield(pln, 'propDoseCalc') && isfield(pln.propDoseCalc,'useWaterPhantom')
                    this.useWaterPhantom = pln.propDoseCalc.useWaterPhantom;
                else
                    this.useWaterPhantom = false;
                end

                if isfield(pln, 'propDoseCalc') && isfield(pln.propDoseCalc,'roomMaterial')
                    this.roomMaterial = pln.propDoseCalc.roomMaterial;
                else
                    this.roomMaterial = 'Air';
                end
                 
             end

            this.FREDrootFolder = fullfile(matRad_cfg.matRadRoot, 'FRED');

            tmpWlsFolder = this.FREDrootFolder;
            tmpWlsFolder(this.FREDrootFolder == '\') = '/';
            this.fragDataLib = [tmpWlsFolder, '/CarbonFragmentationLibraries/data'];

        end

        function set.FREDrootFolder(obj, pathValue)
            obj.FREDrootFolder = pathValue;

            obj.updatePaths;
        end
    end

    methods(Access = protected)

        dij = calcDose(this,ct,cst,stf)


        function dij = initDoseCalc(this,ct,cst,stf)

            matRad_cfg = MatRad_Config.instance();

            dij = initDoseCalc@DoseEngines.matRad_MonteCarloEngineAbstract(this,ct,cst,stf); 
            
            dij = this.allocateQuantityMatrixContainers(dij, {'physicalDose'});

            
            %!!!!!!!!!! This is a mess !!!!! TODO: Understand this better
            if any(strcmp(this.machine.meta.machine, 'generic_MCsquare'))
%                this.machine.meta.fitWithSpotSizeAirCorrection = false;
                this.machine.meta.fitWithSpotSizeAirCorrection = false;

            else
                this.machine.meta.fitWithSpotSizeAirCorrection = true;
            end

            %%% just for testing
            if isempty(this.numHistoriesDirect)
                this.numHistoriesDirect = matRad_cfg.propDoseCalc.defaultNumHistoriesDirect;
            end


            if isempty(this.numHistoriesPerBeamlet)
                this.numHistoriesPerBeamlet = matRad_cfg.propDoseCalc.defaultNumHistoriesPerBeamlet;
            end

            
            %Issue a warning when we have more than 1 scenario
            if dij.numOfScenarios ~= 1
                matRad_cfg.dispWarning('FRED is only implemented for single scenario use at the moment. Will only use the first Scenario for Monte Carlo calculation!');
            end
            
            % prefill ordering of MCsquare bixels
            %dij.FREDCalcOrder = NaN*ones(dij.totalNumOfBixels,1);  
            
            % Check for model consistency
            if ~isempty(this.bioParam)
                %try to load the machine
                % if isa(this.bioParam, 'matRad_BiologicalModel')
                %     this.bioParam.machine = this.machine;
                % end

%                this.calcBioDose = this.bioParam.calcBioDose;
                this.calcBioDose = this.bioParam.bioOpt;

                % if any(strcmp(this.bioParam.RequiredBaseData, {'LET'}))
                %     this.calcLET = true;
                % end
            else
                this.calcBioDose = 0;
            end
            
            if this.calcBioDose

                dij = this.loadBiologicalBaseData(cst,dij);
                dij = this.allocateQuantityMatrixContainers(dij,{'mAlphaDose', 'mSqrtBetaDose'});
            end

            % if isa(this.bioParam, 'matRad_bioModel_constRBE')
            %     dij.RBE = this.bioParam.RBE;
            % end
            if strcmp(this.bioParam.model, 'constRBE')
                dij.RBE = this.bioParam.RBE;
            end
            
            
            
            if this.calcLET
                if this.calcDoseDirect
                    this.scorers = [this.scorers, {'LETd'}];

                    dij = this.allocateQuantityMatrixContainers(dij, {'mLETDose'});
                else
                    this.calcLET = false;

                    matRad_cfg.dispWarning('LETd Dij calculation not yet implemented');
                end
                %matRad_cfg.dispWarning('LET calculation not yet supported');
            end

            %if isempty(this.useInternalHUConversion)
            %    this.useInternalHUConversion = false;
                %For the time being and for testing just use the internal
                %one
                %this.useInternalHUConversion = true;
            %end

        end

       function writeTreeDirectory(this)

            %fred_cfg = MatRad_FREDConfig.instance();

            if ~exist(this.MCrunFolder, 'dir')
                mkdir(this.MCrunFolder);
            end

            %write input folder
            if ~exist(this.inputFolder, 'dir')
                mkdir(this.inputFolder);
            end

            %Build MCrun/inp/regions and
            %      MCrun/inp/plan
            if ~exist(this.regionsFolder, 'dir')
                mkdir(this.regionsFolder);
            end

            if ~exist(this.planFolder, 'dir')
                mkdir(this.planFolder);
            end
        end

        writeRunFile(~, fName)

        %% Write files functions
                
        writeRegionsFile(this,fName, stf)

        writePlanDeliveryFile(this, fName, stf)
  
        writePlanFile(this,fName, stf)

        function writeFredInputAllFiles(this,stf)
    
            %fred_cfg = MatRad_FREDConfig.instance();
            
            %write fred.inp file
            runFilename = fullfile(this.MCrunFolder, this.runInputFilename);
            this.writeRunFile(runFilename);
        
            %write region/region.inp file
            regionFilename = fullfile(this.regionsFolder, this.regionsFilename);
            this.writeRegionsFile(regionFilename, stf);
        
            %write plan file
            planFile = fullfile(this.planFolder, this.planFilename);
            this.writePlanFile(planFile,stf);

            %write planDelivery file

            planDeliveryFile = fullfile(this.planFolder,this.planDeliveryFilename);
            this.writePlanDeliveryFile(planDeliveryFile, stf);
        end

        function dij = loadBiologicalBaseData(this,cst,dij)
            matRad_cfg = MatRad_Config.instance();
            
            matRad_cfg.dispInfo('Initializing biological dose calculation...\n');
            
            dij.ax              = zeros(dij.doseGrid.numOfVoxels,1);
            dij.bx              = zeros(dij.doseGrid.numOfVoxels,1);
            
            cstDownsampled = matRad_setOverlapPriorities(cst);
            
            % resizing cst to dose cube resolution
            cstDownsampled = matRad_resizeCstToGrid(cstDownsampled,dij.ctGrid.x,dij.ctGrid.y,dij.ctGrid.z,...
                dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);
            % retrieve photon LQM parameter for the current dose grid voxels
            [dij.ax,dij.bx] = matRad_getPhotonLQMParameters(cstDownsampled,dij.doseGrid.numOfVoxels,1,this.VdoseGrid);
            
            this.vAlphaX = dij.ax;
            this.vBetaX = dij.bx;
        end


        function dij = allocateQuantityMatrixContainers(this,dij,names)

            if this.calcDoseDirect
                numOfBixelsContainer = 1;
            else
                numOfBixelsContainer = dij.totalNumOfBixels;
            end
            
            %Loop over all requested quantities
            for n = 1:numel(names)
                %Create Cell arrays for container and dij
                szContainer = [numOfBixelsContainer size(this.multScen.scenMask)];
                %tmpMatrixContainers.(names{n}) = cell(szContainer);
                dij.(names{n}) = cell(size(this.multScen.scenMask));
                
                %Now preallocate a matrix in each active scenario using the
                %scenmask
                if this.calcDoseDirect
                    dij.(names{n})(this.multScen.scenMask) = {zeros(dij.doseGrid.numOfVoxels,this.numOfColumnsDij)};
                else
                    %We preallocate a sparse matrix with sparsity of
                    %1e-3 to make the filling slightly faster
                    %TODO: the preallocation could probably
                    %have more accurate estimates
                    dij.(names{n})(this.multScen.scenMask) = {spalloc(dij.doseGrid.numOfVoxels,this.numOfColumnsDij,round(prod(dij.doseGrid.numOfVoxels,this.numOfColumnsDij)*1e-3))};
                end
            end

        end


    end

    methods (Static)

        function [available,msg] = isAvailable(pln,machine)   
            % see superclass for information
            
            msg = [];
            available = false;

            if nargin < 2
                machine = matRad_loadMachine(pln);
            end

            %checkBasic
            try
                checkBasic = isfield(machine,'meta') && isfield(machine,'data');

                %check modality
                checkModality = any(strcmp(DoseEngines.matRad_ParticleFREDEngine.possibleRadiationModes, machine.meta.radiationMode));
                
                preCheck = checkBasic && checkModality;

                if ~preCheck
                    return;
                end
            catch
                msg = 'Your machine file is invalid and does not contain the basic field (meta/data/radiationMode)!';
                return;
            end


            %For the time being, just use the generic machine, will need to
            %have a specific one later on
            switch machine.meta.radiationMode
                case 'protons'
                    available = any(strcmp(pln.machine,{'Generic', 'generic_MCsquare', 'newGeneric_4Aug', 'test_machine_SAD'}));

                case 'carbon'
                    available = any(strcmp(pln.machine,{'Generic', 'HITfixedBL'}));
            end
            msg = 'Machine check is currently not reliable';

        end

        function cube = readSimulationOutput(this,fileName,doseGrid)

            if exist('doseGrid', 'var')
                resolution = doseGrid;
            end
            cube = matRad_readMhd(fileName);
        end

        
        function dijMatrix = readSparseDijBinOlderVersion(fName)
 
            f = fopen(fName,'r','l');

            %Header
            dims = fread(f,3,"int32");
            res = fread(f,3,"float32");
            offset = fread(f,3,"float32");
            numberOfBixels = fread(f,1,"int32");
        
            values = [];
            voxelIndices = [];
            colIndices = [];
            
            fprintf("Reading %d number of beamlets in %d voxels (%dx%dx%d)\n",numberOfBixels,prod(dims),dims(1),dims(2),dims(3));
        
            for i = 1:numberOfBixels
                %Read Beamlet
                bixNum = fread(f,1,"int32");
                beamNum = fread(f,1,"int32");
                numVox  = fread(f,1,"int32");
                
                colIndices(end+1:end+numVox) = bixNum + 1;
                currVoxelIndices = fread(f,numVox,"uint32") + 1;
                values(end+1:end+numVox)  = fread(f,numVox,"float32");
            
                [indX, indY, indZ] = ind2sub(dims, currVoxelIndices);

                voxelIndices(end+1:end+numVox) = sub2ind(dims([2,1,3]), indY, indX, indZ);
                fprintf("\tRead beamlet %d, Field %d, %d voxels...\n",bixNum,beamNum,numVox);
            end
            
            fclose(f);
            
            dijMatrix = sparse(voxelIndices,colIndices,values,prod(dims),numberOfBixels);
        end

        function dijMatrix = readSparseDijBin(fName)

            f = fopen(fName,'r','l');

            %Header
            dims = fread(f,3,"int32");
            %dims = dims([1,2,3]);
            res = fread(f,3,"float32");
            offset = fread(f,3,"float32");
            numberOfBixels = fread(f,1,"int32");
        
            values = [];
            voxelIndices = [];
            colIndices = [];
            
            fprintf("Reading %d number of beamlets in %d voxels (%dx%dx%d)\n",numberOfBixels,prod(dims),dims(1),dims(2),dims(3));
        
            for i = 1:numberOfBixels
                %Read Beamlet
                bixNum = fread(f,1,"int32");
                numVox  = fread(f,1,"int32");
                
                colIndices(end+1:end+numVox) = bixNum + 1;
                currVoxelIndices = fread(f,numVox,"uint32") + 1;
                values(end+1:end+numVox)  = fread(f,numVox,"float32");
            
                [indX, indY, indZ] = ind2sub(dims, currVoxelIndices);

%                voxelIndices(end+1:end+numVox) = sub2ind(dims, indY, indX, indZ);
                voxelIndices(end+1:end+numVox) = sub2ind(dims([2,1,3]), indY, indX, indZ);

                fprintf("\tRead beamlet %d, %d voxels...\n",bixNum,numVox);
            end
            
            fclose(f);
            
            dijMatrix = sparse(voxelIndices,colIndices,values,prod(dims),numberOfBixels);

        end


    end


     methods (Access = private)

        function updatePaths(obj)
            obj.MCrunFolder     = fullfile(obj.FREDrootFolder, 'MCrun');
            obj.inputFolder     = fullfile(obj.MCrunFolder, 'inp');
            obj.regionsFolder   = fullfile(obj.inputFolder, 'regions');
            obj.planFolder      = fullfile(obj.inputFolder, 'plan');
        end


        function calculationAvailable = checkSystemAvailability(this)

            calculationAvailable = true;

 
            matRad_cfg = MatRad_Config.instance();

            % Disabling GPU for WSL
            if this.useWSL
                matRad_cfg.dispInfo('Calling FRED from WSL');

                if this.useGPU
                    matRad_cfg.dispWarning('GPU not available in WSL');
                end
                this.useGPU = false;
            end

            
            switch this.machine.meta.radiationMode
                case 'protons'

                case 'carbon'
                    % This is temporary
                    if ~this.useWSL
                        matRad_cfg.dispWarning('carbon simulation only available in WSL mode. Switching to WSL');
                        this.useWSL = true;
                    end
                    
                    availableVersionsForCarbon = {'3.69.14'};

                    if ~ismember(this.currentVersion, availableVersionsForCarbon)
                        matRad_cfg.dispWarning('current version: %s does not support carbon calculation', this.currentVersion);
                        calculationAvailable = false;
                    else
                        %this only works in wsl
                        this.createSymbolicLinkToData();
                    end

            end

            this.availableVersions = this.getAvailableVersions;

        end


        function version = getVersion(this)

            matRad_cfg = MatRad_Config.instance();

            [status, cmdOut] = system([this.cmdCall,'fred -vn']);

            if status == 0
                %dotIdx = find('.' == cmdOut, 1,'first')-1;
                version = cmdOut(1:end-1);

            else
                matRad_dispError('Something wrong occured in checking FRED installation. Please check correct FRED installation');
            end
        end

        function availableVersions = getAvailableVersions(this, sist)
            
            
            matRad_cfg = MatRad_Config.instance();

            if ~exist('sist', 'var') || isempty(sist)
                if this.useWSL
                    sist = 'wsl';
                else
                    sist = 'win';
                end                   
            end

            if strcmp(sist, 'win')
                currCmdCall = '';
            else
                currCmdCall = 'wsl if [ -f ~/.fredenv.sh ] ; then source ~/.fredenv.sh ; fi;';
            end

            availableVersions = [];

            [status, cmdOut] = system([currCmdCall, 'fred -listVers']);
            
            if status == 0
                nLidx = regexp(cmdOut, '\n')+6; %6 because of tab
                nVersions = numel(nLidx)-1;
                
                for versIdx=1:nVersions
                    availableVersions = [availableVersions,{cmdOut(nLidx(versIdx):nLidx(versIdx)+5)}];
                end

            else
                matRad_cfg.dispError('Something wrong occured in checking FRED available verions. Please check correct FRED installation');
            end
        end
        
        function createSymbolicLinkToData(this)
            matRad_cfg = MatRad_Config.instance();
            % This only works for wsl
            originDataDirectory = this.fragDataLib;
            originDataDirectory(1:find(originDataDirectory == ':')) = [];
            originDataDirectory = ['/mnt/c', originDataDirectory];

            tmpWlsFolder = this.MCrunFolder;
            tmpWlsFolder(this.MCrunFolder == '\') = '/';

            targetDataDirectory = tmpWlsFolder;
            targetDataDirectory(1:find(targetDataDirectory == ':')) = [];
            targetDataDirectory = ['/mnt/c', targetDataDirectory];

            currCmdCall = sprintf('ln -s %s %s',originDataDirectory,targetDataDirectory);

            
            
            
            
            
            
            if matRad_cfg.logLevel > 1
                [stat,~] = system([this.cmdCall, currCmdCall], '-echo');
            else    
                [stat,~] = system([this.cmdCall, currCmdCall]);
            end
            
            
            
            if stat > 0
                rmCmdCode = sprintf('rm %s/data', targetDataDirectory);
                if matRad_cfg.logLevel>1
                    [~,~] = system(['wsl ', rmCmdCode], '-echo');
                    [~,~] = system([this.cmdCall, currCmdCall], '-echo');
                else
                    [~,~] = system(['wsl ', rmCmdCode]);
                    [~,~] = system([this.cmdCall, currCmdCall]);

                end
            end
        end

     
     end

     

     methods

         function set.sourceModel(this, value)
            matRad_cfg = MatRad_Config.instance();

            valid = ischar(value) && any(strcmp(value, this.AvailableSourceModels));

            if valid
                this.sourceModel = value;
            end
               
         end

         function version = get.currentVersion(this)

             if isempty(this.currentVersion)
                version = this.getVersion();
                this.currentVersion = version;
             else
                version = this.currentVersion;
             end
               
         end

         function set.useWSL(this,value)
             
             this.useWSL = value;
        
             if value
                this.cmdCall = 'wsl if [ -f ~/.fredenv.sh ] ; then source ~/.fredenv.sh ; fi;';
             else
                this.cmdCall = '';
             end
         end

     end
end

