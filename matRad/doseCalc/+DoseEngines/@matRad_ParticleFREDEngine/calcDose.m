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
    
    currFolder = pwd;

    cd(this.FREDrootFolder);

    %Now we can run initDoseCalc as usual
    dij = this.initDoseCalc(ct,cst,stf);

    HUcube{1} =  matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y',  dij.ctGrid.z,ct.cubeHU{1}, ...
                                dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'linear');

     switch this.HUtable
        case 'internal'
            if any(HUcube{1}(:)>this.hLutLimits(2)) || any(HUcube{1}(:)<this.hLutLimits(1))
                matRad_cfg.dispWarning('HU outside of boundaries');
                this.HUclamping = true;
            end
        otherwise
            matRad_cfg.dispInfo('Using custom HU table: %s\n', this.HUtable);
     end
     
    %Write the directory tree necessary for the simulation
    this.writeTreeDirectory();
    
    % Write patient
    cd(this.regionsFolder);
    
    %%%%%% !!!!!!!!!!!!! mind this flip !!!!!!!!!!!!! %%%%%
    HUcube{1} = flip(HUcube{1},2);
    
    fileNamePatient = fullfile(this.regionsFolder, this.patientFilename);
    patientMetadata.imageOrigin = [0 0 0];
    patientMetadata.resolution  = [this.doseGrid.resolution.x, this.doseGrid.resolution.y, this.doseGrid.resolution.z];
    patientMetadata.datatype = 'int16';
    matRad_writeMHD(fileNamePatient, HUcube{1},patientMetadata);

    cd(this.FREDrootFolder);

    getPointAtBAMS = @(target,source,distance,BAMStoIso) (target -source)*(-BAMStoIso)/distance + source;%(target  + source*(BAMStoIso - distance))/distance;
          
    % Loop over the stf to rearrange data
    counter = 0;

    % Attention here the field in amchine.meta has been artificially set in
    % the initDoseCalc
    emittanceBaseData = matRad_MCemittanceBaseData(this.machine,stf);

    if this.calcDoseDirect
        cumulativeWeights = 0;
        for i=1:length(stf)
            cumulativeWeights = cumulativeWeights + sum([stf(i).ray.weight]);
        end
    
    end
    
    for i = 1:length(stf)

        stfFred(i).gantryAngle     = stf(i).gantryAngle;
        stfFred(i).couchAngle      = stf(i).couchAngle;
        
        % Get dose grid resolution
        doseGridResolution = [this.doseGrid.resolution.x, this.doseGrid.resolution.y, this.doseGrid.resolution.z];
        
        % get isocenter coordinates in dose cube coordinate system.
        isoInDoseGridCoord = matRad_world2cubeCoords(stf(i).isoCenter,this.doseGrid);

        % First voxel in dose cube system has center at resolution. Thus
        % surface is at 1/2 resolution wrt zero of that system
        fredCubeSurfaceInDoseCubeCoords =  0.5*doseGridResolution;
        
        % Get coordinates of pivot point in FRED cube (center of
        % geometrical cube) in dose cube coordinates.
        fredPivotInCubeCoordinates = 0.5*this.doseGrid.dimensions.*doseGridResolution + fredCubeSurfaceInDoseCubeCoords;
        
        % Place isocenter wrt to FRED coordinate system. FRED then:
        % - place component pivot in te center of FRED coordinate system
        % - translate the component by ISO coordinartes, so that ISO point
        %   is at the center of FRED world coordinates.
        % - apply gantry and couch rotations
        stfFred(i).isoCenter = -(fredPivotInCubeCoordinates - isoInDoseGridCoord);

        % First coordinate is flipped
        stfFred(i).isoCenter = stfFred(i).isoCenter.*[-1 1 1];
        

        nominalEnergies        = unique([stf(i).ray.energy]);
        [~,nominalEnergiesIdx] = intersect([this.machine.data.energy],nominalEnergies);
        
        energyIdxInEmittance   = ismember(emittanceBaseData.energyIndex, nominalEnergiesIdx);
        monteCarloBaseData     = emittanceBaseData.monteCarloData(energyIdxInEmittance);
        

        stfFred(i).nominalEnergies = nominalEnergies;
        stfFred(i).energies        = [monteCarloBaseData.MeanEnergy].*this.numOfNucleons./this.primaryMass;
        stfFred(i).energySpread    = [monteCarloBaseData.EnergySpread];
        stfFred(i).energySpreadMeV = [monteCarloBaseData.EnergySpread].*[monteCarloBaseData.MeanEnergy]/100;
        stfFred(i).FWHMs = 2.355*[monteCarloBaseData.SpotSize1x];
        
        stfFred(i).energySpreadFWHMMev    = 2.355*stfFred(i).energySpreadMeV;
        stfFred(i).BAMStoIsoDist   = emittanceBaseData.nozzleToIso;
        
        % This is just used for reference and define the field size
        %stfFred(i).FWHMs = 2.355*[monteCarloBaseData.SpotSize1x(1)];
        switch this.sourceModel

            case 'gaussian'
        
            case 'emittance'
                stfFred(i).emittanceX = [];
                stfFred(i).twissBetaX  = [];
                stfFred(i).twissAlphaX = [];
                stfFred(i).emittanceRefPlaneDistance = [];
                
                % Need to get the parameters for the model from MCemittance
                for eIdx=emittanceBaseData.energyIndex
                    % Only using first focus index for now
                    tmpOpticsData = emittanceBaseData.fitBeamOpticsForEnergy(eIdx,1);
                    stfFred(i).emittanceX       = [stfFred(i).emittanceX,  tmpOpticsData.twissEpsilonX];
                    stfFred(i).twissBetaX       = [stfFred(i).twissBetaX,  tmpOpticsData.twissBetaX];
                    stfFred(i).twissAlphaX      = [stfFred(i).twissAlphaX, tmpOpticsData.twissAlphaX];
                    stfFred(i).emittanceRefPlaneDistance = [stfFred(i).emittanceRefPlaneDistance, this.machine.meta.BAMStoIsoDist];
                end

            case 'sigmaSqrModel'
                stfFred(i).sSQr_a = [];
                stfFred(i).sSQr_b = [];
                stfFred(i).sSQr_c = [];

                for eIdx=emittanceBaseData.energyIndex
                    tmpOpticsData = emittanceBaseData.fitBeamOpticsForEnergy(eIdx,1);
                    stfFred(i).sSQr_a           = [stfFred(i).sSQr_a, tmpOpticsData.sSQ_a];
                    stfFred(i).sSQr_b           = [stfFred(i).sSQr_b, tmpOpticsData.sSQ_b];
                    stfFred(i).sSQr_c           = [stfFred(i).sSQr_c, tmpOpticsData.sSQ_c];
                end
            otherwise
                matRad_cfg.dispWarning('Unrecognized source model, setting gaussian');

        end

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
                dij.beamNum(counter,1)  = i;
                dij.rayNum(counter,1)   = j;
                dij.bixelNum(counter,1) = k;
            end
           
            for k = 1:numel(stfFred(i).energies)

                if any(stf(i).ray(j).energy == stfFred(i).nominalEnergies(k)) %any(stf(i).ray(j).energy == nominalEnergies(k))
                    stfFred(i).energyLayer(k).rayNum   = [stfFred(i).energyLayer(k).rayNum j];
                    
                    if isfield(stf(i).ray(j), 'energySpread')
                    %     % If also energy spread is provided, need to
                    %     % match both energy and energy spread. This is
                    %     % a special case, different energy spreads are
                    %     % seen as different energy layers
                        stfFred(i).energyLayer(k).bixelNum = [stfFred(i).energyLayer(k).bixelNum ...
                                                              find(stf(i).ray(j).energy == stfFred(i).nominalEnergies(k) & ...
                                                              stf(i).ray(j).energySpread == stfFred(i).energySpread(k))];

                    else

                        stfFred(i).energyLayer(k).bixelNum = [stfFred(i).energyLayer(k).bixelNum ...
                            find(stf(i).ray(j).energy == stfFred(i).nominalEnergies(k))];
                    end
                    targetX = stf(i).ray(j).targetPoint_bev(1);
                    targetY = stf(i).ray(j).targetPoint_bev(3);

                    sourceX = stf(i).ray(j).rayPos_bev(1);
                    sourceY = stf(i).ray(j).rayPos_bev(3);

                    distance = stf(i).ray(j).targetPoint_bev(2) - stf(i).ray(j).rayPos_bev(2);
                    
                    divergenceX = (targetX - sourceX)/distance;
                    divergenceY = (targetY - sourceY)/distance;
                    
                    stfFred(i).energyLayer(k).targetPoints    = [stfFred(i).energyLayer(k).targetPoints; -targetX targetY];

                    stfFred(i).energyLayer(k).rayPosX         = [stfFred(i).energyLayer(k).rayPosX, getPointAtBAMS(targetX,sourceX,distance,stfFred(i).BAMStoIsoDist)];
                    stfFred(i).energyLayer(k).rayPosY         = [stfFred(i).energyLayer(k).rayPosY, getPointAtBAMS(targetY,sourceY,distance,stfFred(i).BAMStoIsoDist)];

                    stfFred(i).energyLayer(k).rayDivX         = [stfFred(i).energyLayer(k).rayDivX, divergenceX];
                    stfFred(i).energyLayer(k).rayDivY         = [stfFred(i).energyLayer(k).rayDivY, divergenceY];
                    
                    
                    if this.calcDoseDirect
                            stfFred(i).energyLayer(k).numOfPrimaries = [stfFred(i).energyLayer(k).numOfPrimaries ...
                                                  stf(i).ray(j).weight(stf(i).ray(j).energy == stfFred(i).nominalEnergies(k))];
                    else

                         stfFred(i).energyLayer(k).numOfPrimaries = [stfFred(i).energyLayer(k).numOfPrimaries, 1];
                    end

                end

            end
        end

        %FRED works in cm
        stfFred(i).isoCenter       = stfFred(i).isoCenter/10;
        stfFred(i).BAMStoIsoDist   = stfFred(i).BAMStoIsoDist/10;

        switch this.sourceModel
            case 'gaussian'
                stfFred(i).FWHMs           = stfFred(i).FWHMs/10;                
            case 'emittance'
                stfFred(i).emittanceRefPlaneDistance = stfFred(i).emittanceRefPlaneDistance/10;
            case 'sigmaSqrModel'

        end

        stfFred(i).totalNumOfBixels = stf(i).totalNumOfBixels;
        for j=1:numel(stfFred(i).nominalEnergies)
           stfFred(i).energyLayer(j).rayPosX      = stfFred(i).energyLayer(j).rayPosX/10;
           stfFred(i).energyLayer(j).rayPosY      = stfFred(i).energyLayer(j).rayPosY/10;
           stfFred(i).energyLayer(j).targetPoints = stfFred(i).energyLayer(j).targetPoints/10;
           stfFred(i).energyLayer(j).nBixels      = numel(stfFred(i).energyLayer(j).bixelNum);
           %This is necessary because of the sum(w) at the end of
           %calcDoseForward
           if this.calcDoseDirect
               stfFred(i).energyLayer(j).numOfPrimaries = this.conversionFactor*stfFred(i).energyLayer(j).numOfPrimaries/cumulativeWeights;
           end
        end
    end
    
    counterFred = 0;
    FredOrder = NaN * ones(dij.totalNumOfBixels,1);
    for i = 1:length(stf)
        for j = 1:numel(stfFred(i).nominalEnergies) %numel(stfFred(i).energies)
            for k = 1:numel(stfFred(i).energyLayer(j).numOfPrimaries)
                counterFred = counterFred + 1;
                ix = find(i                               == dij.beamNum & ...
                    stfFred(i).energyLayer(j).rayNum(k)   == dij.rayNum & ...
                    stfFred(i).energyLayer(j).bixelNum(k) == dij.bixelNum);
    
                FredOrder(ix) = counterFred;
            end
        end
    end
    
    if any(isnan(FredOrder))
        matRad_cfg.dispError('Invalid ordering of Beamlets for FRED computation!');
    end
    
    % %% MC computation and dij filling

    this.writeFredInputAllFiles(stfFred);

    if ~this.exportCalculation
        %should add checks for installation of fred and so on

        % Check consistency of installation
        if this.checkSystemAvailability()
            
            matRad_cfg.dispInfo('calling FRED');

            cd(this.MCrunFolder);

            if ~this.useGPU
                cmdCall = [this.cmdCall, 'fred -f fred.inp -nogpu'];
            else
                cmdCall = [this.cmdCall, 'fred -f fred.inp'];
            end

            % printOutput to matLab console
            if this.printOutput
                [status,~] = system(cmdCall,'-echo');
            else
                [status,~] = system(cmdCall);
            end
            cd(this.FREDrootFolder);
        else
            matRad_cfg.dispError('FRED setup incorrect for this plan simulation');
        end

        if status==0
            matRad_cfg.dispInfo(' done\n');
        end


       if ~this.calcDoseDirect

            doseDijFolder = fullfile(this.MCrunFolder, 'out', 'scoreij');
            doseDijFile = 'Phantom.Dose.bin';
            dijFileName = fullfile(doseDijFolder,doseDijFile);
            dijMatrix = this.readSparseDijBin(dijFileName);

            % Check consistency
            if isequal(size(dijMatrix), [dij.doseGrid.numOfVoxels,dij.totalNumOfBixels])
                %When scoring dij, FRED internaly normalizes to 1
                dij.physicalDose{1} = this.conversionFactor*dijMatrix(:,FredOrder);
            end

            if this.calcLET
                dijFile = 'Phantom.LETd.bin';
                dijFileName = fullfile(doseDijFolder,dijFile);
                dijMatrix = this.readSparseDijBin(dijFileName);
                
                if isequal(size(dijMatrix), [dij.doseGrid.numOfVoxels,dij.totalNumOfBixels])
                    % Need to divide by 10, FRED scores in /cm (?)
                    dij.mLETd{1} = dijMatrix(:,FredOrder)./10;
                end
                
                dij.mLETDose{1} = sparse(dij.physicalDose{1}.*dij.mLETd{1});

            end

        else

            %if str2num(this.currentVersion(1:find(this.currentVersion == '.',1,'last')-1)) > 3.6
            doseCubeFolder = fullfile(this.MCrunFolder, 'out', 'score');
            doseCubeFile = 'Phantom.Dose.mhd';
            %else
            %         doseCubeFolder = fullfile(this.MCrunFolder, 'out');
            %         doseCubeFile = 'Dose.mhd';
            %end
            
            cube = matRad_readMhd(doseCubeFolder, doseCubeFile);

            % readMHD internaly flips dimension 2 of the cube. IDK why this
            % is done, probably required by MCsquare. But with FRED's
            % coordinate system there is no need for such a flip, so for the
            % time being we'll just revert it here
            cube = flip(cube,2);

            dij.physicalDose{1} = sparse(this.VdoseGrid, ones(numel(this.VdoseGrid),1),cube(this.VdoseGrid), dij.doseGrid.numOfVoxels,1);

            if this.calcLET
                cubeLET = matRad_readMhd(fullfile(this.MCrunFolder, 'out', 'score'), 'Phantom.LETd.mhd');
                
                % same as cube
                cubeLET = flip(cubeLET,2);

                dij.mLETd{1}    = sparse(this.VdoseGrid, ones(numel(this.VdoseGrid),1),cubeLET(this.VdoseGrid)./10, dij.doseGrid.numOfVoxels,1);
                dij.mLETDose{1} = sparse(this.VdoseGrid, ones(numel(this.VdoseGrid),1),(cubeLET(this.VdoseGrid)./10).*cube(this.VdoseGrid), dij.doseGrid.numOfVoxels,1);

            end

       end

       if this.calcBioDose
           % recover alpha and beta maps
           % radDepths not used by models
           tmpBixel.radDepths = zeros(size(this.VdoseGrid,1),1);
        
           tmpBixel.vAlphaX   = dij.ax{1}(this.VdoseGrid);
           tmpBixel.vBetaX    = dij.bx{1}(this.VdoseGrid);
           tmpBixel.vABratio  = dij.ax{1}(this.VdoseGrid)./dij.bx{1}(this.VdoseGrid);

           if this.calcDoseDirect
                tmpKernel.LET = dij.mLETd{1}(this.VdoseGrid);

                [tmpBixelAlpha, tmpBixelBeta] = this.bioParam.calcLQParameterForKernel(tmpBixel,tmpKernel);
                
                tmpBixelAlpha(isnan(tmpBixelAlpha)) = 0;
                tmpBixelBeta(isnan(tmpBixelBeta)) =  0;

                dij.mAlphaDose{1}     = sparse(this.VdoseGrid, ones(numel(this.VdoseGrid),1),tmpBixelAlpha.*dij.physicalDose{1}(this.VdoseGrid), dij.doseGrid.numOfVoxels,1);
                dij.mSqrtBetaDose{1} = sparse(this.VdoseGrid, ones(numel(this.VdoseGrid),1),sqrt(tmpBixelBeta).*dij.physicalDose{1}(this.VdoseGrid), dij.doseGrid.numOfVoxels,1);
           else    
               % Loop over all bixels
               for bxlIdx = 1:dij.totalNumOfBixels
                   bixelLET           = full(dij.mLETd{1}(:,bxlIdx));
                   tmpKernel.LET      = bixelLET(this.VdoseGrid);

                   [tmpBixelAlpha, tmpBixelBeta] = this.bioParam.calcLQParameterForKernel(tmpBixel,tmpKernel);

                   tmpBixelAlpha(isnan(tmpBixelAlpha)) = 0;
                   tmpBixelBeta(isnan(tmpBixelBeta)) =  0;
                
                   dij.mAlphaDose{1}(:,bxlIdx)     = sparse(this.VdoseGrid, ones(numel(this.VdoseGrid),1),tmpBixelAlpha.*dij.physicalDose{1}(this.VdoseGrid,bxlIdx), dij.doseGrid.numOfVoxels,1);
                   dij.mSqrtBetaDose{1}(:,bxlIdx) = sparse(this.VdoseGrid, ones(numel(this.VdoseGrid),1),sqrt(tmpBixelBeta).*dij.physicalDose{1}(this.VdoseGrid,bxlIdx), dij.doseGrid.numOfVoxels,1);
               end
           end
        end 
    else
        matRad_cfg.dispInfo('All files have been generated');
    end
    
    dij = this.finalizeDose(dij);
    

    cd(currFolder);

end