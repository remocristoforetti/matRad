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

        constantRBE = NaN;              % constant RBE value
        useInternalHUConversion;
        noozleToAxis;
        scorers = {'Dose'};
        calcLET;
     
        HUclamping = false;
        HUtable;
        defaultHUtable = 'matRad_water.inp';
        defaultHUtable = 'matRad_water.txt';
    end

    properties
        exportCalculation = false;
    properties (SetAccess = private)
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
            
            if nargin < 1
                pln = [];
            end

            % call superclass constructor
            this = this@DoseEngines.matRad_MonteCarloEngineAbstract(pln);

            % check if bio optimization is needed and set the
            % coresponding boolean accordingly
            % TODO:
            % This should not be handled here as an optimization property
            % We should rather make optimization dependent on what we have
            % decided to calculate here.
            % if nargin > 0 
            %     if (isfield(pln,'propOpt')&& isfield(pln.propOpt,'bioOptimization')&& ...
            %         (isequal(pln.propOpt.bioOptimization,'LEMIV_effect') ||...
            %         isequal(pln.propOpt.bioOptimization,'LEMIV_RBExD')) && ...
            %         strcmp(pln.radiationMode,'carbon'))
            %     this.calcBioDose = true;
            %     elseif strcmp(pln.radiationMode,'protons') && isfield(pln,'propOpt') && isfield(pln.propOpt,'bioOptimization') && isequal(pln.propOpt.bioOptimization,'const_RBExD')
            %         this.constantRBE = 1.1;                    
            %     end
            % end
            if nargin > 0 
                if (isfield(pln,'propOpt')&& isfield(pln.propOpt,'bioOptimization')&& ...
                    (isequal(pln.propOpt.bioOptimization,'LEMIV_effect') ||...
                    isequal(pln.propOpt.bioOptimization,'LEMIV_RBExD')) && ...
                    strcmp(pln.radiationMode,'carbon'))
                this.calcBioDose = flase;
                matRad_cfg.dispWarning('bio calculation not yet implemented in this version.');
                elseif strcmp(pln.radiationMode,'protons') && isfield(pln,'propOpt') && isfield(pln.propOpt,'bioOptimization') && isequal(pln.propOpt.bioOptimization,'const_RBExD')
                    this.constantRBE = NaN;
                     matRad_cfg.dispWarning('bio calculation not yet implemented in this version.');
                end
            end

            %Instantiate a MatRad_FREDConfig
            %This instance only contains properties and default values that
            %can be changed at any time
            
            matRad_cfg = MatRad_Config.instance();
            fred_cfg = MatRad_FREDConfig.instance();
            
            fred_cfg.FREDrootFolder = fullfile(matRad_cfg.matRadRoot, 'FRED');
        end                
            this.FREDrootFolder = fullfile(matRad_cfg.matRadRoot, 'FRED');
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
            
            %%% just for testing
            %this.numHistoriesDirect = 1000000;
            %this.numHistoriesPerBeamlet = 100000;


            %Issue a warning when we have more than 1 scenario
            if dij.numOfScenarios ~= 1
                matRad_cfg.dispWarning('FRED is only implemented for single scenario use at the moment. Will only use the first Scenario for Monte Carlo calculation!');
            end
            
            % prefill ordering of MCsquare bixels
            dij.FREDCalcOrder = NaN*ones(dij.totalNumOfBixels,1);  
            
            if ~isnan(this.constantRBE) 
                dij.RBE = this.constantRBE;
            end
            

            if ~this.calcDoseDirect
                this.calcDoseDirect = true;
                matRad_cfg.dispWarning('Dij calculation still not supported. Setting to directDoseCalc')
            end
            
            if this.calcLET 
                %this.scorers = [this.scorers, {'LETd'}];
                this.calcLET = false;
                matRad_cfg.dispWarning('LET calculation not yet supported');
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
        
        
        % function printArray(~,fID,arrayName, array, arrayElementType)
        % 
        %     fprintf(fID, arrayName);
        %     for k=1:numel(array)-1
        %         fprintf(fID, [arrayElementType, ','], array(k));
        %     end
        %     fprintf(fID, [arrayElementType, ']'], array(end));
        % end

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
            available = any(strcmp(pln.machine,{'Generic', 'generic_MCsquare'}));
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
end

