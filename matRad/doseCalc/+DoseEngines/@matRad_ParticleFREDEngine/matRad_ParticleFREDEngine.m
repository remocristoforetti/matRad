classdef matRad_ParticleFREDEngine < DoseEngines.matRad_MonteCarloEngineAbstract
% Engine for particle dose calculation using monte carlo calculation
% specificly the FRED MC code ()
% for more informations see superclass
% DoseEngines.matRad_MonteCarloEngineAbstract
%
% pln.propDoseCalc fields:
% 
% HUclamping:            [b] allows for clamping of HU table
% HUtable:               [s] 'internal', 'matRad_water'
% exportCalculation      [b] t: Only write simulation paramter files
%                            f: run FRED
% sourceModel            [s] see AvailableSourceModels, {'gaussian', 'emittance', 'sigmaSqrModel'}
% useWSL                 [b]
% useGPU                 [b]
% roomMaterial           [s] vacuum, Air
% printOutput            [b]
% numHistoriesDirect     [i]
% numHistoriesPerBeamlet [i]


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
        
        possibleRadiationModes = {'protons'};
        name = 'FRED';
        shortName = 'FRED';
    end

    properties (SetAccess = public, GetAccess = public) %(SetAccess = protected, GetAccess = public)
        
        currFolder = pwd; %folder path when set

        noozleToAxis;
        scorers;
     
        HUtable;
        defaultHUtable = 'matRad_default_FRED';
        AvailableSourceModels = {'gaussian', 'emittance', 'sigmaSqrModel'};
        calcBioDose;

        sourceModel;
        

        useWaterPhantom;
        waterPhantomIpot;
        roomMaterial;
        useWSL;
        useGPU;
        currentVersion;
        availableVersions = {'3.70.0'};
        fragDataLib;
        printOutput;
        radiationMode;
        primaryMass;
        numOfNucleons;
    end

    properties
        exportCalculation = false;
        calcLET;
        constantRBE = NaN;              % constant RBE value
        HUclamping;

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
            %   pln: matRad plan meta information struct
            
            matRad_cfg = MatRad_Config.instance();
            if nargin < 1
                pln = [];
            end

            % call superclass constructor
            this = this@DoseEngines.matRad_MonteCarloEngineAbstract(pln);

            if nargin > 0
                if isfield(pln, 'radiationMode')
                    this.radiationMode = pln.radiationMode;
                end
            end

            this.FREDrootFolder = fullfile(matRad_cfg.primaryUserFolder, 'FRED');
            
            if ~exist(this.FREDrootFolder, 'dir')
                mkdir(this.FREDrootFolder);
            end

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

            %Issue a warning when we have more than 1 scenario
            if dij.numOfScenarios ~= 1
                matRad_cfg.dispWarning('FRED is only implemented for single scenario use at the moment. Will only use the first Scenario for Monte Carlo calculation!');
            end
            

            % Check for model consistency
            if ~isempty(this.bioParam)
                this.calcBioDose = this.bioParam.bioOpt;
            else
                this.calcBioDose = 0;
            end
            
            % Limit RBE calculation to proton models for the time being
            if this.calcBioDose

                switch this.radiationMode

                    case 'protons'

                        dij = this.loadBiologicalBaseData(cst,dij);
                        dij = this.allocateQuantityMatrixContainers(dij,{'mAlphaDose', 'mSqrtBetaDose'});

                        % Only considering LET based models 
                        %dij.calcLET   = true;
                        this.calcLET = true;
                    otherwise
                        matRad_cfg.dispWarning('biological dose calculation not supported for radiation modality: %s', this.radiationMode);
                        this.calcBioDose = false;
                end
            end

            if strcmp(this.bioParam.model, 'constRBE')
                dij.RBE = this.bioParam.RBE;
            end
            
            if this.calcLET
                this.scorers = [this.scorers, {'LETd'}];
                % Allocate containers for both LET*Dose and dose weighted
                % LET. This last is used for biological calculation as well
                dij = this.allocateQuantityMatrixContainers(dij, {'mLETDose', 'mLETd'});
            end

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
            [dij.ax,dij.bx] = matRad_getPhotonLQMParameters(cstDownsampled,dij.doseGrid.numOfVoxels,this.VdoseGrid);
            
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

    methods

        function setDefaults(this)
            setDefaults@DoseEngines.matRad_MonteCarloEngineAbstract(this);

            matRad_cfg = MatRad_Config.instance();

            this.HUclamping             = false;
            this.HUtable                = this.defaultHUtable;
            this.exportCalculation      = false;
            this.sourceModel            = this.AvailableSourceModels{1};
            this.useWSL                 = false;
            this.useGPU                 = true;
            this.roomMaterial           = 'Air';
            this.printOutput            = true;
            this.numHistoriesDirect     = matRad_cfg.defaults.propDoseCalc.numHistoriesDirect;
            this.numHistoriesPerBeamlet = matRad_cfg.defaults.propDoseCalc.numHistoriesPerBeamlet;
            this.scorers                = {'Dose'};
            this.primaryMass            = 1.00727;
            this.numOfNucleons          = 1;



        end

        function writeHlut(this,hLutFile)

            matRad_cfg = MatRad_Config.instance();
            fileName = fullfile(this.regionsFolder, 'hLut.inp');

            mainFolder        = fullfile(matRad_cfg.matRadSrcRoot,'hluts');
            userDefinedFolder = fullfile(matRad_cfg.primaryUserFolder, 'hluts');

            % Collect all the subfolders
            if ispc
                folderDelimiter = ';';
            elseif isunix
                folderDelimiter = ':';
            else
                folderDelimiter = ';';
            end

            searchPath = [strsplit(genpath(mainFolder), folderDelimiter)';strsplit(genpath(userDefinedFolder), folderDelimiter)'];
            searchPath(cellfun(@isempty, searchPath)) = [];
            
            availableHLUTs = cellfun(@(x) dir([x,'\*.hlut']), searchPath, 'UniformOutput',false);
            availableHLUTs = cell2mat(availableHLUTs);

            hLUTindex = find(strcmp([hLutFile,'.hlut'], {availableHLUTs.name}));

            if ~isempty(hLUTindex)
                selectedHlutfile = fullfile(availableHLUTs(hLUTindex).folder, availableHLUTs(hLUTindex).name);
                
                template = fileread(selectedHlutfile);
        
                newLut = fopen(fileName, 'w');
                fprintf(newLut, template);
                fclose(newLut);
            else
                matRad_cfg.dispError(['Cannot open hLut: ',hLutFile]);
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

            available = preCheck;
        end

        function cube = readSimulationOutput(this,folder,fileName,doseGrid)

            if exist('doseGrid', 'var')
                resolution = doseGrid;
            end
            cube = matRad_readMhd(folder,fileName);
        end

        function dijMatrix = readSparseDijBin(fName)

            f = fopen(fName,'r','l');

            %Header
            fileFormatVerison = fread(f,1,"int32");
            dims = fread(f,3,"int32");
            res = fread(f,3,"float32");
            offset = fread(f,3,"float32");
            nComponents = fread(f,1,"int32");
            numberOfBixels = fread(f,1,"int32");
        
            values = [];
            valuesDen = [];
            voxelIndices = [];
            colIndices = [];
            valuesNom = [];
            
            fprintf("Reading %d number of beamlets in %d voxels (%dx%dx%d)\n",numberOfBixels,prod(dims),dims(1),dims(2),dims(3));
        
            for i = 1:numberOfBixels
                %Read Beamlet
                bixNum = fread(f,1,"int32");
                numVox  = fread(f,1,"int32");
                
                colIndices(end+1:end+numVox) = bixNum + 1;
                currVoxelIndices = fread(f,numVox,"uint32") + 1;
                tmpValues = fread(f,numVox*nComponents,"float32");
                valuesNom = tmpValues(1:nComponents:end);%tmpValues(nComponents:nComponents:end);
                % values(end+1:end+numVox) = tmpValues(1:nComponents:end);%tmpValues(nComponents:nComponents:end);

                if nComponents == 2
                    valuesDen = tmpValues(nComponents:nComponents:end);
                    values(end+1:end+numVox) = valuesNom./valuesDen;
                else
                    values(end+1:end+numVox) = valuesNom;
                end
            
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

                otherwise
                    matRad_cfg.dispError('Only proton dose calculation available with this version of FRED');
                    calculationAvailable = false;

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

        function [radiationMode] = updateRadiationMode(this,value)
            matRad_cfg = MatRad_Config.instance();
            
            if any(strcmp(value, this.possibleRadiationModes))
                radiationMode = value;
            else
                matRad_cfg.dispError('Invalid radiation modality: %s', value);
            end

            switch radiationMode
                case 'protons'
                    this.primaryMass = 1.00727; % Da
                    this.numOfNucleons = 1;
                otherwise
                    matRad_cfg.dispError('Only proton dose calculation available with this version of FRED');

            end
            
            matRad_cfg.dispWarning('Selected radiation modality: %s with primary mass: %2.3f', radiationMode, this.primaryMass);
        end
     
     end

     

     methods

         function set.sourceModel(this, value)
            matRad_cfg = MatRad_Config.instance();

            valid = ischar(value) && any(strcmp(value, this.AvailableSourceModels));

            if valid
                
                this.sourceModel = value;
            else
                matRad_cfg.dispWarning('Unable to set source model:%s, setting default:%s', value, this.AvailableSourceModels{1})
                this.sourceModel = this.AvailableSourceModels{1};
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

         function set.radiationMode(this,value)
            if ischar(value)
                if ~isempty(this.radiationMode) && ~strcmp(this.radiationMode, value)
                    this.radiationMode = value;
                    this.updateRadiationMode(this.radiationMode);
                elseif isempty(this.radiationMode)
                    this.radiationMode = value;
                end
            end

         end

     end
end

