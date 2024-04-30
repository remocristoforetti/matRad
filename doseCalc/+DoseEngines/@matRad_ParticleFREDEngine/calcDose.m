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
    %fred_cfg = MatRad_FREDConfig.instance();
    %matRad_cfg.dispError('Still need to implement this function');
    
    currFolder = pwd;

    % cd to FRED folder (necessary ?)

    cd(this.FREDrootFolder);

    %Now we can run initDoseCalc as usual
    dij = this.initDoseCalc(ct,cst,stf);

     for s = 1:dij.numOfScenarios
         HUcube{s} =  matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y',  dij.ctGrid.z,ct.cubeHU{s}, ...
                                     dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'linear');
     end

     if this.useInternalHUConversion
        if any(HUcube{1}(:)>this.hLutLimits(2)) || any(HUcube{1}(:)<this.hLutLimits(1))
            matRad_cfg.dispWarning('HU outside of boundaries');
            this.HUclamping = true;
        end
        this.HUtable = 'internal';
     else
        this.HUtable = this.defaultHUtable;
     end
   

    %Write the directory tree necessary for the simulation
    this.writeTreeDirectory();
    

    if ~this.useWaterPhantom
        cd(this.regionsFolder);
        matRad_writeMhd(HUcube{1},[this.doseGrid.resolution.x, this.doseGrid.resolution.y, this.doseGrid.resolution.z], this.patientFilename, 'MET_SHORT');
        cd(this.FREDrootFolder);
    end

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
        
        % stfFred(i).isoCenter       = [-0.5*dij.doseGrid.resolution.x*(dij.doseGrid.dimensions(2))+stf(i).isoCenter(1),...
        %                               0.5*dij.doseGrid.resolution.y*(dij.doseGrid.dimensions(1))-stf(i).isoCenter(2),...
        %                               0.5*dij.doseGrid.resolution.z*(dij.doseGrid.dimensions(3))-stf(i).isoCenter(3)];
        % 
        % 
        % isoOffset = [-dij.doseGrid.resolution.x dij.doseGrid.resolution.y dij.doseGrid.resolution.z]./2;
        % 
        % stfFred(i).isoCenter = (-stfFred(i).isoCenter + isoOffset - [-ct.resolution.x, ct.resolution.y, ct.resolution.z]);

        stfFred(i).isoCenter       = [0.5*dij.doseGrid.resolution.x*(dij.doseGrid.dimensions(2))-stf(i).isoCenter(1),...
                                      0.5*dij.doseGrid.resolution.y*(dij.doseGrid.dimensions(1))-stf(i).isoCenter(2),...
                                      0.5*dij.doseGrid.resolution.z*(dij.doseGrid.dimensions(3))-stf(i).isoCenter(3)];

        
        isoOffset = [dij.doseGrid.resolution.x dij.doseGrid.resolution.y dij.doseGrid.resolution.z]./2;

        stfFred(i).isoCenter = (-stfFred(i).isoCenter + isoOffset - [ct.resolution.x, ct.resolution.y, ct.resolution.z]);
        stfFred(i).isoCenter(1) = -stfFred(i).isoCenter(1);

        nominalEnergies        = unique([stf(i).ray.energy]);
        [~,nominalEnergiesIdx] = intersect([this.machine.data.energy],nominalEnergies);
        
        energyIdxInEmittance   = ismember(emittanceBaseData.energyIndex, nominalEnergiesIdx);

        monteCarloBaseData     = emittanceBaseData.monteCarloData(energyIdxInEmittance);
        
        % TODO: Find a better setup, this is only used to test multiuple Es
        % for baseData analysis
        if isfield(stf(i).ray, 'energySpread')
            %Only works if only one ray with multiple energies is defined
            nEnergies = numel(stf(i).ray.energy);
            stfFred(i).nominalEnergies = stf(i).ray.energy;
            stfFred(i).energies = ones(1,nEnergies)*monteCarloBaseData.MeanEnergy;
            stfFred(i).energySpread = stf.ray.energySpread;
            stfFred(i).energySpreadMeV = (stf.ray.energySpread*monteCarloBaseData.MeanEnergy/100);        
            stfFred(i).FWHMs = 2.355*[monteCarloBaseData.SpotSize1x(stf(i).ray(1).focusIx(1))]*ones(1,numel(stfFred(i).energies));

        else
            stfFred(i).nominalEnergies = nominalEnergies;
            stfFred(i).energies        = [monteCarloBaseData.MeanEnergy];
            stfFred(i).energySpread    = [monteCarloBaseData.EnergySpread];
            stfFred(i).energySpreadMeV = [monteCarloBaseData.EnergySpread].*[monteCarloBaseData.MeanEnergy]/100;
            stfFred(i).FWHMs = 2.355*[monteCarloBaseData.SpotSize1x(stf(i).ray(1).focusIx(1))];

        end
        
        
        % This energy spread is sigmaE (?) this should match MC2 better
        %stfFred(i).energySpread    = 2.355*energySpreadMeV + 0.1;
        stfFred(i).energySpreadFWHMMev    = 2.355*stfFred(i).energySpreadMeV;
        stfFred(i).BAMStoIsoDist   = emittanceBaseData.nozzleToIso;%this.machine.meta.BAMStoIsoDist;
        
        % This is just used for reference and define the field size
        %stfFred(i).FWHMs = 2.355*[monteCarloBaseData.SpotSize1x(1)];
        switch this.sourceModel

            case 'gaussian'
        
            case 'emittance'
                stfFred(i).emittanceX       = [monteCarloBaseData.twissEpsilonX];
                stfFred(i).twissBetaX       = [monteCarloBaseData.twissBetaX];
                stfFred(i).twissAlphaX      = [monteCarloBaseData.twissAlphaX];

            case 'sigmaSqrModel'
                stfFred(i).sSQr_a           = [monteCarloBaseData.sSQ_a];
                stfFred(i).sSQr_b           = [monteCarloBaseData.sSQ_b];
                stfFred(i).sSQr_c           = [monteCarloBaseData.sSQ_c];
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
                dij.beamNum(counter)  = i;
                dij.rayNum(counter)   = j;
                dij.bixelNum(counter) = k;
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
                    
                    %Normalization not needed in principle, this is
                    %handled internally by FRED

                    % normDivergence = 1;% (divergenceX + divergenceY);
                    % divergenceX = divergenceX/normDivergence;
                    % divergenceY = divergenceY/normDivergence;



                    stfFred(i).energyLayer(k).targetPoints    = [stfFred(i).energyLayer(k).targetPoints; -targetX targetY];

                    stfFred(i).energyLayer(k).rayPosX         = [stfFred(i).energyLayer(k).rayPosX, getPointAtBAMS(targetX,sourceX,distance,stfFred(i).BAMStoIsoDist)];
                    stfFred(i).energyLayer(k).rayPosY         = [stfFred(i).energyLayer(k).rayPosY, getPointAtBAMS(targetY,sourceY,distance,stfFred(i).BAMStoIsoDist)];

                    stfFred(i).energyLayer(k).rayDivX         = [stfFred(i).energyLayer(k).rayDivX, divergenceX];
                    stfFred(i).energyLayer(k).rayDivY         = [stfFred(i).energyLayer(k).rayDivY, divergenceY];
                    
                    
                    if this.calcDoseDirect

                        %     stfFred(i).energyLayer(k).numOfPrimaries = [stfFred(i).energyLayer(k).numOfPrimaries ...
                        %                           stf(i).ray(j).weight(stf(i).ray(j).energy == stfFred(i).nominalEnergies(k) & ...
                        %                           stf(i).ray(j).energySpread == stfFred(i).energySpread(k))];
                        % 
                        % else
                            stfFred(i).energyLayer(k).numOfPrimaries = [stfFred(i).energyLayer(k).numOfPrimaries ...
                                                  stf(i).ray(j).weight(stf(i).ray(j).energy == stfFred(i).nominalEnergies(k))];
                        % end

                    else
                         %matRad_cfg.dispWarning('DIJ calculation not yet implemented');
                         stfFred(i).energyLayer(k).numOfPrimaries = [stfFred(i).energyLayer(k).numOfPrimaries ...
                             1];
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
                %stfFred(i).refPlane          = stfFred(i).refPlane/10;
                %stfFred(i).emittanceX       = stfFred(i).emittanceX./10;
                %stfFred(i).twissBetaX       = stfFred(i).twissBetaX./10;
            case 'sigmaSqrModel'

        end

        stfFred(i).totalNumOfBixels = stf(i).totalNumOfBixels;
        for j=1:numel(stfFred(i).nominalEnergies) %size(stfFred(i).energies,2)
           stfFred(i).energyLayer(j).rayPosX      = stfFred(i).energyLayer(j).rayPosX/10;
           stfFred(i).energyLayer(j).rayPosY      = stfFred(i).energyLayer(j).rayPosY/10;
           stfFred(i).energyLayer(j).targetPoints = stfFred(i).energyLayer(j).targetPoints/10;
           stfFred(i).energyLayer(j).nBixels      = numel(stfFred(i).energyLayer(j).bixelNum); %numel(stfFred(i).energyLayer(j).rayPosX);
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
                ix = find(i                                   == dij.beamNum & ...
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
        matRad_cfg.dispInfo('calling FRED');


        %Need to make this better
        if this.checkSystemAvailability()
            cd(this.MCrunFolder);
            if matRad_cfg.logLevel>1

                [status,~] = system([this.cmdCall, 'fred -f fred.inp'],'-echo');
            else
                [status,~] = system([this.cmdCall, 'fred -f fred.inp']);
            end
            cd(this.FREDrootFolder);
        else
            matRad_cfg.dispError('FRED setup incorrect for this plan simulation');
        end

        if status==0
            matRad_cfg.dispInfo('done\n');
        end

        if ~this.calcDoseDirect

            %Need to add sanity check for presence of Dij.bin
            %Now working for one field. Need to check what happens with two

            switch this.currentVersion
                case '3.69.14'
                    doseDijFolder = fullfile(this.MCrunFolder, 'out', 'scoreij');
                    doseDijFile = 'Phantom.Dose.bin';
                    dijFileName = fullfile(doseDijFolder,doseDijFile);
                    dijMatrix = this.readSparseDijBin(dijFileName);

                otherwise
                    doseDijFolder = fullfile(this.MCrunFolder, 'out');
                    doseDijFile = 'Dij.bin';
                    dijFileName = fullfile(doseDijFolder,doseDijFile);
                    dijMatrix = this.readSparseDijBinOlderVersion(dijFileName);

            end

            % dijFileName = fullfile(doseDijFolder,doseDijFile);
            % 
            % %This introduces a permutation in the collected cube
            % dijMatrix = this.readSparseDijBin(dijFileName);

            % Check consistency
            if isequal(size(dijMatrix), [dij.doseGrid.numOfVoxels,dij.totalNumOfBixels])
                %When scoring dij, FRED internaly normalizes to 1
                dij.physicalDose{1} = this.conversionFactor*dijMatrix(:,FredOrder);
            end

        else
            % For the time being, just load the dose and resample

            switch this.currentVersion
                case '3.69.14'
                    doseCubeFolder = fullfile(this.MCrunFolder, 'out', 'score');
                    doseCubeFile = 'Phantom.Dose.mhd';
                otherwise
                    doseCubeFolder = fullfile(this.MCrunFolder, 'out');
                    doseCubeFile = 'Dose.mhd';
            end

            cube = matRad_readMhd(doseCubeFolder, doseCubeFile);

            % readMHD internaly flips dimension 2 of the cube. IDK why this
            % is done, probably required by MCsquare. But with FRED's
            % coordinate system there is no need of sich a flip, so for the
            % time being we'll just revert it here
            cube = flip(cube,2);

            dij.physicalDose{1} = sparse(this.VdoseGrid, ones(numel(this.VdoseGrid),1),cube(this.VdoseGrid), dij.doseGrid.numOfVoxels,1);

            if this.calcLET
                cubeLET = matRad_readMhd(fullfile(this.MCrunFolder, 'out'), 'LETd.mhd');
                
                % same as cube
                cubeLET = flip(cubeLET,2);

                dij.mLETDose{1} = sparse(this.VdoseGrid, ones(numel(this.VdoseGrid),1),cubeLET(this.VdoseGrid).*cube(this.VdoseGrid), dij.doseGrid.numOfVoxels,1);
            end

            
            if this.calcBioDose
                 
                % recover alpha and beta maps
                tmpBixel.radDepths = zeros(dij.doseGrid.numOfVoxels,1);
                tmpBixel.vAlphaX   = this.vAlphaX;
                tmpBixel.vBetaX    = this.vBetaX;
                tmpBixel.vABratio   = this.vAlphaX./this.vBetaX;
                tmpKernel.LET       =  cubeLET(:);


                [tmpBixelAlpha, tmpBixelBeta] = this.bioParam.calcLQParameterForKernel(tmpBixel,tmpKernel);

                tmpBixelAlpha(isnan(tmpBixelAlpha)) = 0;
                tmpBixelBeta(isnan(tmpBixelBeta)) =  0;
                dij.mAlphaDose = {tmpBixelAlpha.*dij.physicalDose{1}};
                dij.mSqrtBetaDose = {sqrt(tmpBixelBeta).*dij.physicalDose{1}};

                %bioDoseCube = matRad_readMhd(fullfile(this.MCrunFolder, 'out', 'RBE'), ['DoseBio_', fName, '.mhd']);
                %dij.BioDose{1} = sparse(this.VdoseGrid, ones(numel(this.VdoseGrid),1),bioDoseCube(this.VdoseGrid), dij.doseGrid.numOfVoxels,1);
            end


        end
    else
        matRad_cfg.dispInfo('All files have been generated');
    end
    
    dij = this.finalizeDose(dij);
    

    cd(currFolder);

end