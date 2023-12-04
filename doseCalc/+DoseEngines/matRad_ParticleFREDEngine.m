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
        phantomPatientName = 'CTpatient.mhd'
        scorers = {'Dose'};
        calcLET;
        runSimulationFileName = 'fred.inp'
        noozleToAxis;
        exportCalculation = true;
        %regionsFolder = '/inp/regions'
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
        end                
    end
    
    methods(Access = protected)
        
        function dij = calcDose(this,ct,cst,stf)
            % matRad FRED monte carlo proton dose calculation wrapper
            % can be automaticly called through matRad_calcDose or
            % matRad_calcParticleDoseMC
            %
            %
            % call
            %   dij = this.calcDose(ct,stf,pln,cst)
            %
            % input
            %   ct:          	matRad ct struct
            %   cst:            matRad cst struct
            %   stf:         	atRad steering information struct
            %   
            % output
            %   dij:            matRad dij struct
            %
            % References
            %
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
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            matRad_cfg = MatRad_Config.instance();

            %matRad_cfg.dispError('Still need to implement this function');
            
            this.currFolder = pwd;
            fullfilename = mfilename('fullpath');
            FredFolder = [matRad_cfg.matRadRoot filesep 'FRED'];

            % cd to FRED folder (necessary ?)
            cd(FredFolder);

            %Now we can run initDoseCalc as usual
            dij = this.initDoseCalc(ct,cst,stf);

            %check for presence of HUcorrection file
            if ~exist([FredFolder filesep 'HUmaterialConversionTables'],'dir')
                this.useInternalHUConversion = true;
            end

            %For the time being, use FRED with the dose grid resolution, will deal
            %with this later
            for s = 1:dij.numOfScenarios
                HUcube{s} =  matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y',  dij.ctGrid.z,ct.cubeHU{s}, ...
                                            dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'linear');
            end

            %Write the directory tree necessary for teh simulation
            this.writeTreeDirectory();

            %write mhd patient cube
            %patientFilename = [pwd,filesep, 'MCrun\inp\regions', filesep, this.phantomPatientName];
            currTmpFolder = pwd;
            cd([FredFolder, filesep, 'MCrun', filesep, 'inp', filesep, 'regions\']);
            matRad_writeMhd(HUcube{1},[this.doseGrid.resolution.x, this.doseGrid.resolution.y, this.doseGrid.resolution.z], this.phantomPatientName, 'MET_SHORT');
            cd(currTmpFolder);

            getPointAtBAMS = @(target,source,distance,BAMStoIso) (target  + source*(BAMStoIso - distance))/distance;
            % Loop over the stf to rearrange data
            counter = 0;

            for i = 1:length(stf)

                stfFred(i).gantryAngle = stf(i).gantryAngle;
                stfFred(i).couchAngle      = stf(i).couchAngle;
                stfFred(i).isoCenter       = stf(i).isoCenter;
                stfFred(i).energies        = unique([stf(i).ray.energy]);
                stfFred(i).BAMStoIsoDist   = this.machine.meta.BAMStoIsoDist;

                [~,eIdx] = intersect(stfFred(i).energies, [this.machine.data.energy]);

                %for the time being use just initial sigma, then need to
                %go through MCemittanceBaseData
                stfFred(i).FWHMs = [];
                for j= eIdx'
                    stfFred(i).FWHMs           = [stfFred(i).FWHMs, [this.machine.data(j).initFocus.sigma(1)]];
                end
                stfFred(i).FWHMs = 2.355*stfFred(i).FWHMs;

                % allocate empty target point container
                for j = 1:numel(stfFred(i).energies)
                    stfFred(i).energyLayer(j).targetPoints   = [];
                    stfFred(i).energyLayer(j).numOfPrimaries = [];
                    stfFred(i).energyLayer(j).rayNum         = [];
                    stfFred(i).energyLayer(j).bixelNum       = [];
                    stfFred(i).energyLayer(j).rayDivX        = [];
                    stfFred(i).energyLayer(j).rayDivY        = [];
                    stfFred(i).energyLayer(j).rayPosX        = [];
                    stfFred(i).energyLayer(j).rayPosY        = [];
                end

                for j = 1:stf(i).numOfRays
                    for k = 1:stf(i).numOfBixelsPerRay(j)
                        counter = counter + 1;
                        dij.beamNum(counter)  = i;
                        dij.rayNum(counter)   = j;
                        dij.bixelNum(counter) = k;
                    end
             
                    
                    for k = 1:numel(stfFred(i).energies)

                        if any(stf(i).ray(j).energy == stfFred(i).energies(k))
                            stfFred(i).energyLayer(k).rayNum   = [stfFred(i).energyLayer(k).rayNum j];
                            stfFred(i).energyLayer(k).bixelNum = [stfFred(i).energyLayer(k).bixelNum ...
                                find(stf(i).ray(j).energy == stfFred(i).energies(k))];

                            targetX = -stf(i).ray(j).targetPoint_bev(1);
                            targetY = stf(i).ray(j).targetPoint_bev(3);

                            sourceX = -stf(i).ray(j).rayPos_bev(1);
                            sourceY = stf(i).ray(j).rayPos_bev(3);

                            distance = stf(i).ray(j).targetPoint_bev(2) - stf(i).ray(j).rayPos_bev(2);
                            
                            divergenceX = (targetX - sourceX)/distance;
                            divergenceY = (targetY - sourceY)/distance;
                            
                            %Normalization not needed in principle, this is
                            %handled internally by FRED

                            normDivergence = 1;% (divergenceX + divergenceY);
                            divergenceX = divergenceX/normDivergence;
                            divergenceY = divergenceY/normDivergence;


                            stfFred(i).energyLayer(k).targetPoints    = [stfFred(i).energyLayer(k).targetPoints; targetX targetY];
                            %stfFred(i).energyLayer(k).targetPointDist = [stfFred(i).energyLayer(k).targetPointDist; distance];

                            stfFred(i).energyLayer(k).rayPosX         = [stfFred(i).energyLayer(k).rayPosX, getPointAtBAMS(targetX,sourceX,distance,stfFred(i).BAMStoIsoDist)];
                            stfFred(i).energyLayer(k).rayPosY         = [stfFred(i).energyLayer(k).rayPosY, getPointAtBAMS(targetY,sourceY,distance,stfFred(i).BAMStoIsoDist)];

                            stfFred(i).energyLayer(k).rayDivX         = [stfFred(i).energyLayer(k).rayDivX, divergenceX];
                            stfFred(i).energyLayer(k).rayDivY         = [stfFred(i).energyLayer(k).rayDivY, divergenceY];
                            
                            if this.calcDoseDirect
                                 stfFred(i).energyLayer(k).numOfPrimaries = [stfFred(i).energyLayer(k).numOfPrimaries ...
                                                      stf(i).ray(j).weight(stf(i).ray(j).energy == stfFred(i).energies(k))];
                            else
                                 matRad_cfg.dispWarning('DIJ calculation not yet implemented');
                                 stfFred(i).energyLayer(k).numOfPrimaries = [stfFred(i).energyLayer(k).numOfPrimaries ...
                                     MCsquareConfig.Num_Primaries];
                            end

                        end
                        
                    end
                end
 
                stfFred(i).nBixels = [];
 
                for k=1:numel(stfFred(i).energies)
                    stfFred(i).nBixels = [stfFred(i).nBixels, numel(stfFred(i).energyLayer(k).rayNum)];
                end
                %FRED works in cm
                stfFred(i).isoCenter       = stfFred(i).isoCenter/10;
                stfFred(i).BAMStoIsoDist   = stfFred(i).BAMStoIsoDist/10;
                stfFred(i).FWHMs           = stfFred(i).FWHMs/10;
                 
                for j=1:size(stfFred(i).energies,2)
                   stfFred(i).energyLayer(j).rayPosX      = stfFred(i).energyLayer(j).rayPosX/10;
                   stfFred(i).energyLayer(j).rayPosY      = stfFred(i).energyLayer(j).rayPosY/10;
                   stfFred(i).energyLayer(j).targetPoints = stfFred(i).energyLayer(j).targetPoints/10;
                end
            end

             % remember order

             counterFred = 0;
             FredOrder = NaN * ones(dij.totalNumOfBixels,1);
             for i = 1:length(stf)
                 for j = 1:numel(stfFred(i).energies)
                     for k = 1:numel(stfFred(i).energyLayer(j).numOfPrimaries)
                         counterFred = counterFred + 1;
                         ix = find(i                                     == dij.beamNum & ...
                                   stfFred(i).energyLayer(j).rayNum(k)   == dij.rayNum & ...
                                   stfFred(i).energyLayer(j).bixelNum(k) == dij.bixelNum);

                         FredOrder(ix) = counterFred;
                     end
                 end
             end
             
            if any(isnan(FredOrder))
                matRad_cfg.dispError('Invalid ordering of Beamlets for MCsquare computation!');
            end
             
            % %% MC computation and dij filling

            this.writeFredInputAllFiles(stfFred);

            if ~this.exportCalculation
                matRad_cfg.dispError('Auto calculation not yet implemented');
            else
                matRad_cfg.dispInfo('All files have been generated');
                dij.physicalDose = [];
            end
            
            this.finalizeDose(ct,cst,stf,dij);
            
            cd(this.currFolder);
        end

        function dij = initDoseCalc(this,ct,cst,stf)
            matRad_cfg = MatRad_Config.instance();            
            
            dij = initDoseCalc@DoseEngines.matRad_MonteCarloEngineAbstract(this,ct,cst,stf); 
            
            %%% just for testing
            this.numHistoriesDirect = 10000;
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
                matRad_cfg.dispWarning('Dij calculation still not supported. Setting ')
            end
            
            if this.calcLET 
                %this.scorers = [this.scorers, {'LETd'}];
                this.calcLET = false;
                matRad_cfg.dispWarning('LET calculation not yet supported');
            end

        end

       function writeTreeDirectory(this)

            %FRED forder already set as current directory
            runFolder = [pwd, filesep, 'MCrun'];

            if ~exist(runFolder, 'dir')

                mkdir(runFolder);
            end

            % Build input folder, /MCrun/inp
            inputFolder = [runFolder, filesep, 'inp'];
            
            if ~exist(inputFolder, 'dir')
                mkdir(inputFolder);
            end

            %Build MCrun/inp/regions and
            %      MCrun/inp/plan
            regionsFolder = [inputFolder, filesep, 'regions'];
            planFolder    = [inputFolder, filesep, 'plan'];

            if ~exist(regionsFolder, 'dir')
                mkdir(regionsFolder);
            end
            if ~exist(planFolder, 'dir')
                mkdir(planFolder);
            end
        end

         function writeRunFile(~, fName)
            matRad_cfg = MatRad_Config.instance();

            fID = fopen(fName, 'w');

            try
                fprintf(fID, 'include: inp/funcs.inp\n');
                fprintf(fID, 'include: inp/regions/regions.inp\n');
                fprintf(fID, 'include: inp/plan/planDelivery.inp\n');
            catch
                matRad_cfg.dispError('Failed to write run file');
            end

            fclose(fID);
        end
        
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

        function writeFredInputAllFiles(this,stf)
            
            %write fred.inp file
            runFilename = [pwd, filesep, 'MCrun', filesep, this.runSimulationFileName];
            this.writeRunFile(runFilename);

            %write region/region.inp file
            regionFilename = [pwd, filesep, 'MCrun', filesep, 'inp', filesep, 'regions', filesep, 'regions.inp'];
            this.writeRegionsFile(regionFilename);
            
            %write funcs, this is here untill we figure out how to solve
            %this mess
            funcsFile = [pwd, filesep, 'MCrun', filesep,'inp', filesep, 'funcs.inp'];
            this.writeFuncsFile(funcsFile);

            %write fields file
            fieldsFile = [pwd, filesep, 'MCrun', filesep,'inp', filesep, 'plan', filesep, 'fields.inp'];
            this.writeFieldsFile(fieldsFile, stf);

            %write layers file
            layersFile = [pwd, filesep, 'MCrun', filesep,'inp', filesep, 'plan', filesep, 'layers.inp'];
            this.writeLayersFile(layersFile, stf);

            %write planDelivery file
            planFile = [pwd, filesep, 'MCrun', filesep,'inp', filesep, 'plan', filesep, 'planDelivery.inp'];
            this.writeplanDeliveryFile(planFile, stf);
        end

        function writeRegionsFile(this,fName)
            matRad_cfg = MatRad_Config.instance();

            fID = fopen(fName, 'w');
            
            try
                fprintf(fID,'region<\n');
                fprintf(fID,'\tID=Phantom\n');
                fprintf(fID,'\tCTscan=inp/regions/%s\n', this.phantomPatientName);
                fprintf(fID,'\tpivot=[0.5,0.5,0.5]\n');

                %This is gantry angle = 0
                fprintf(fID, '\tl=[%1.1f,%1.1f,%1.1f]\n', 0,-1,0);
                fprintf(fID, '\tu=[%1.1f,%1.1f,%1.1f]\n', 1,0,0);

                if numel(this.scorers)>1
                    fprintf(fID,'\tscore=[');
                    for k=1:size(this.scorers,2)-1
                        fprintf(fID,'%s,', this.scorers{k});
                    end
                    fprintf(fID,'%s]\n', this.scorers{end});
                else
                    fprintf(fID,'\tscore=[%s]\n', this.scorers{1});
                end

                fprintf(fID,'region>\n');
                if this.useInternalHUConversion
                    fprintf(fID, '\tlUseInternalHU2Mat=t');
                end
            catch
                matRad_cfg.dispError('Failed to write regions file');
            end

            fclose(fID);
        end

        function writeFuncsFile(~, fName)
            matRad_cfg = MatRad_Config.instance();

            fID = fopen(fName, 'w');
            
            try
                fprintf(fID, 'func: valueArr(arr,idx) = arr[idx]\n');
                fprintf(fID, 'func: getDictionaryKeyValue(dic,key) = dic[key]\n');
            catch
                matRad_cfg.dispError('Failed to write funcs file');

            end
            fclose(fID);
        end

        function writeFieldsFile(this,fName,stf)

            matRad_cfg = MatRad_Config.instance();

            fID = fopen(fName, 'w');
            
            % This is hard coded here, could find a cleaner solution


            %Notes
            % The f vector in field is the direction of propagation (by
            % FRED)
            % It is referred to the e3 element of the base, this is for the
            % specific component, thus L(3) is the size of the field in the
            % direction fo propagation. O is referred to the room
            % coordinates system
            try
                fprintf(fID, 'include: inp/plan/layers.inp\n');
                fprintf(fID, 'def: nFields = %i\n',numel(stf));
                for i=1:numel(stf)
                    currFieldID      = i-1;
                    fieldLimX = max(abs([stf(i).energyLayer.rayPosX])) + 10*max([stf(i).FWHMs]);
                    fieldLimY = max(abs([stf(i).energyLayer.rayPosY])) + 10*max([stf(i).FWHMs]);

                    fprintf(fID, 'field: %i; ', currFieldID);
                    fprintf(fID, 'O = [%2.3f,%2.3f,%2.3f]; ',-stf(i).BAMStoIsoDist,0,0);
                    fprintf(fID, 'L = [%2.3f,%2.3f,%2.3f]; ',fieldLimX,fieldLimY,0.1);

                    fprintf(fID, 'pivot = [%0.1f,%0.1f,%0.1f]; ',0.5,0.5,0.5);

                    %Always in z direction, turn the spot and the patient
                    fprintf(fID, 'f = [%i,%i,%i]; ',1,0,0);
                    fprintf(fID, 'u = [%i,%i,%i]\n',0,1,0);
                end
                 fprintf(fID, 'nprim = %i', this.numHistoriesDirect);

            catch
                matRad_cfg.dispError('Failed to write field file');
            end

            fclose(fID);
        end

        function writeLayersFile(this, fName, stf)
            matRad_cfg = MatRad_Config.instance();

            fID = fopen(fName, 'w');
            
            % This is hard coded here, could find a cleaner solution
            try
                layerCounter = 1;
                for i=1:numel(stf)
                    for j=1:size(stf(i).energies,2)

                        fprintf(fID, 'def: L%i = {', layerCounter);
                        fprintf(fID, '"Energy": %3.2f, ', stf(i).energies(j));
                        fprintf(fID, '"FWHM": %3.2f, ',   stf(i).FWHMs(j));
                        
                        %This is the energy spread:
                        fprintf(fID, '"Ts": %f, ',   1);
                        
                        fprintf(fID, '"nSpots": %i, ', stf(i).nBixels(j));
                        
                        % PositionX
                        this.printArray(fID,'"x": [',stf(i).energyLayer(j).rayPosX, '%2.3f');
                        
                        % Position Y
                        this.printArray(fID,', "y": [',stf(i).energyLayer(j).rayPosY, '%2.3f');

                        %divergenceX
                        this.printArray(fID,', "divX": [',stf(i).energyLayer(j).rayDivX, '%2.3f');

                        %divergenceY
                        this.printArray(fID,', "divY": [',stf(i).energyLayer(j).rayDivY, '%2.3f');

                        %w
                        this.printArray(fID,', "w": [',stf(i).energyLayer(j).numOfPrimaries, '%2.3f');

                        fprintf(fID, '}\n');
                        layerCounter = layerCounter+1;
                    end
                end
                
                fprintf(fID, 'def: layers = [');
                for j=1:layerCounter-2
                    fprintf(fID,'L%i,',j);
                end
                    fprintf(fID,'L%i]\n',layerCounter-1);
            catch
                matRad_cfg.dispError('Failed to write layers file');
            end

            fclose(fID);
        end
        

        function writeplanDeliveryFile(this, fName, stf)
            matRad_cfg = MatRad_Config.instance();

            fID = fopen(fName, 'w');
            try
                fprintf(fID, 'include: inp/plan/fields.inp\n');
                fprintf(fID, 'def: ipb = 0\n');
                
                
                %deactivate fileds
                fprintf(fID, 'for(fieldIdx in range(nFields))<\n');
                    fprintf(fID, '\tdeactivate: field_$fieldIdx\n');
                fprintf(fID, 'for>\n\n');

                    %loop over fields
                    fprintf(fID, 'for(fieldIdx in range(nFields))<\n');
                        %activate current filed
                        fprintf(fID, '\tactivate: field_$fieldIdx\n');
    
                        %define pbmaster
                        fprintf(fID, '\tpbmaster: $fieldIdx; Xsec = gauss; columns = [P.x, P.y, N, FWHM, T, v.x, v.y, v.z]\n\n');
                        fprintf(fID, '\tfor(layer in layers)<\n\n');
                        fprintf(fID, '\t\tdef: currEnergy = getDictionaryKeyValue(layer,"Energy")\n');
                        fprintf(fID, '\t\tdef: currFWHM   = getDictionaryKeyValue(layer,"FWHM")\n\n');
                        fprintf(fID, '\t\tfor(pbIdx in range(layer["nSpots"]))<\n');
                            fprintf(fID, '\t\t\tdef: x = valueArr(getDictionaryKeyValue(layer,"x"),pbIdx)\n');
                            fprintf(fID, '\t\t\tdef: y = valueArr(getDictionaryKeyValue(layer,"y"),pbIdx)\n');
                            fprintf(fID, '\t\t\tdef: w = valueArr(getDictionaryKeyValue(layer,"w"),pbIdx)\n');
                            fprintf(fID, '\t\t\tdef: vx = valueArr(getDictionaryKeyValue(layer,"divX"),pbIdx)\n');
                            fprintf(fID, '\t\t\tdef: vy = valueArr(getDictionaryKeyValue(layer,"divY"),pbIdx)\n\n');
    
                            fprintf(fID, '\t\t\tpb: $ipb $fieldIdx $x $y $w $currFWHM $currEnergy $vx $vy 1\n\n');
    
                            fprintf(fID, '\t\t\tdef: ipb = ipb +1;\n');
                        fprintf(fID, '\t\tfor>\n');
                    fprintf(fID, '\tfor>\n');
                    fprintf(fID, '\tdeliver: field_$fieldIdx\n');
                    fprintf(fID, '\tdeactivate: field_$fieldIdx\n');

                    fprintf(fID, 'for>\n\n');
                
            catch
                matRad_cfg.dispError('Failed to write planDelivery file');
            end

            fclose(fID);
        end


        function printArray(~,fID,arrayName, array, arrayElementType)
            
            fprintf(fID, arrayName);
            for k=1:numel(array)-1
                fprintf(fID, [arrayElementType, ','], array(k));
            end
            fprintf(fID, [arrayElementType, ']'], array(end));
        end

    end

    methods (Static)

        % function [binaryFound,binaryFile] = checkBinaries()
        %     % checkBinaries check if the binaries are available on the current
        %     % machine and sets to the mcsquarebinary object property
        %     %        
        %     %
        %     matRad_cfg = MatRad_Config.instance();
        % 
        %     binaryFile = [];
        %     binaryFound = false;
        % 
        %     if ispc
        %         if exist('MCSquare_windows.exe','file') ~= 2
        %             matRad_cfg.dispWarning('Could not find MCsquare binary.\n');
        %         else
        %             binaryFile = 'MCSquare_windows.exe';
        %         end
        %     elseif ismac
        %         if exist('MCsquare_mac','file') ~= 2
        %             matRad_cfg.dispWarning('Could not find MCsquare binary.\n');
        %         else
        %             binaryFile = './MCsquare_mac';
        %         end
        %         %error('MCsquare binaries not available for mac OS.\n');
        %     elseif isunix
        %         if exist('MCsquare_linux','file') ~= 2
        %             matRad_cfg.dispWarning('Could not find MCsquare binary.\n');
        %         else
        %             binaryFile = './MCsquare_linux';
        %         end
        %     end
        % 
        %     if ~isempty(binaryFile)
        %         binaryFound = true;
        %     end
        % 
        % end

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

            %Check the binaries
            %hasBinaries = DoseEngines.matRad_ParticleMCsquareEngine.checkBinaries();
            
            % if ~hasBinaries
            %     return;
            % end

            %For the time being, just use the generic machine, will need to
            %have a specific one later on
            available = strcmp(pln.machine,'Generic');
            msg = 'Machine check is currently not reliable';
        end

        function cube = readSimulationOutput(this,fileName,doseGrid)

            if exist('doseGrid', 'var')
                resolution = doseGrid;
            end
            cube = matRad_readMhd(fileName);
        end
              
    end
end

