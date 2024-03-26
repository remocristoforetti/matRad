classdef matRad_ParticleFREDEngine < DoseEngines.matRad_MonteCarloEngineAbstract
% Engine for particle dose calculation using monte carlo calculation
% specificly the mc square method
% for more informations see superclass
% DoseEngines.matRad_MonteCarloEngineAbstract
%
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

    properties (SetAccess = protected, GetAccess = public)
        
        currFolder = pwd; %folder path when set
        
        %nbThreads; %number of threads for MCsquare, 0 is all available

        useInternalHUConversion;
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
            end

            this.FREDrootFolder = fullfile(matRad_cfg.matRadRoot, 'FRED');

        end

        function set.FREDrootFolder(obj, pathValue)
            obj.FREDrootFolder = pathValue;

            obj.updatePaths;
        end

        % function set.RBEmodel(this, value)
        % 
        %     valid = ischar(value) && any(strcmp(value, this.availableRBEmodels));
        % 
        %     if valid
        %         this.RBEmodel = value;
        %     else
        %         matRad_cfg.dispWarning('RBE model not recognized. Setting constRBE');
        %         this.RBEmodel = 'constRBE';
        %     end
        % 
        % 
        % 
        % end

    end

    methods(Access = protected)

        dij = calcDose(this,ct,cst,stf)


        function dij = initDoseCalc(this,ct,cst,stf)

            matRad_cfg = MatRad_Config.instance();            
            
            dij = initDoseCalc@DoseEngines.matRad_MonteCarloEngineAbstract(this,ct,cst,stf); 
            
            dij = this.allocateQuantityMatrixContainers(dij, {'physicalDose'});

            
            % This is a mess
            if strcmp(this.machine.meta.machine, 'generic_MCsquare')
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

            if isempty(this.useInternalHUConversion)
                this.useInternalHUConversion = false;
                %For the time being and for testing just use the internal
                %one
                %this.useInternalHUConversion = true;
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

        function setUp(this,nCasePerBixel,calcDoseDirect)    
        % SETUP Set up properties used for dose calculation
        %
        % input:
        %   nCasePerBixel:  number of histories per beamlet
        %   calcDoseDirect: binary switch to enable forward dose calculation output
        %
   
            matRad_cfg = MatRad_Config.instance();

            % first argument should be nCasePerBixel
            if (exist('nCasePerBixel','var') && isa(nCasePerBixel,'numeric'))    
                this.nCasePerBixel = nCasePerBixel;
            else
                %set number of particles simulated per pencil beam
                this.nCasePerBixel = matRad_cfg.propMC.MCsquare_defaultHistories;
                matRad_cfg.dispInfo('No number of Histories given or wrong type given. Using default number of Histories per Bixel: %d\n',this.nCasePerBixel);
            end

            if (exist('calcDoseDirect', 'var'))
                this.calcDoseDirect = true;
            end
            
        end
        
        function setBinaries(this)
            % setBinaries check if the binaries are available on the current
            % machine and sets to the mcsquarebinary object property
            %

            [~,binaryFile] = this.checkBinaries();
            this.mcSquareBinary = binaryFile;
        end
        
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
                checkModality = any(strcmp(DoseEngines.matRad_ParticleMCsquareEngine.possibleRadiationModes, machine.meta.radiationMode));
                
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
            available = any(strcmp(pln.machine,{'Generic', 'generic_MCsquare', 'newGeneric_4Aug'}));
            msg = 'Machine check is currently not reliable';
        end

        function cube = readSimulationOutput(this,fileName,doseGrid)

            if exist('doseGrid', 'var')
                resolution = doseGrid;
            end
            cube = matRad_readMhd(fileName);
        end

        
        function dijMatrix = readSparseDijBin(fName)
 
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

                voxelIndices(end+1:end+numVox) = sub2ind(dims, indY, indX, indZ);
                fprintf("\tRead beamlet %d, Field %d, %d voxels...\n",bixNum,beamNum,numVox);
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
     end

     methods

         function set.sourceModel(this, value)
            matRad_cfg = MatRad_Config.instance();

            valid = ischar(value) && any(strcmp(value, this.AvailableSourceModels));

            if valid
                this.sourceModel = value;
            end
               
        end

     end
end

