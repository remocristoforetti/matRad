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
    
    %patientFileName = fullfile(fred_cfg.regionsFolder, fred_cfg.patientFilename);
    cd(this.regionsFolder);
    
    matRad_writeMhd(HUcube{1},[this.doseGrid.resolution.x, this.doseGrid.resolution.y, this.doseGrid.resolution.z], this.patientFilename, 'MET_SHORT');
    cd(this.FREDrootFolder);

    getPointAtBAMS = @(target,source,distance,BAMStoIso) (target -source)*(-BAMStoIso)/distance + source;%(target  + source*(BAMStoIso - distance))/distance;
          
    % Loop over the stf to rearrange data
    counter = 0;

    emittanceBaseData = matRad_MCemittanceBaseData(this.machine,stf);

    cumulativeWeights = 0;
    for i=1:length(stf)
        cumulativeWeights = cumulativeWeights + sum([stf(i).ray.weight]);
    end
    for i = 1:length(stf)

        stfFred(i).gantryAngle     = stf(i).gantryAngle;
        stfFred(i).couchAngle      = stf(i).couchAngle;

        %This is Topas like definition of isocenter
        % stfFred(i).isoCenter       = -[0.5*ct.resolution.x*(ct.cubeDim(2)+1)-stf(i).isoCenter(1),...
        %                               0.5*ct.resolution.y*(ct.cubeDim(1)+1)-stf(i).isoCenter(2),...
        %                               0.5*ct.resolution.z*(ct.cubeDim(3)+1)-stf(i).isoCenter(3)];
        
        stfFred(i).isoCenter       = [0.5*dij.doseGrid.resolution.x*(dij.doseGrid.dimensions(2))-stf(i).isoCenter(1),...
                                      0.5*dij.doseGrid.resolution.y*(dij.doseGrid.dimensions(1))-stf(i).isoCenter(2),...
                                      0.5*dij.doseGrid.resolution.z*(dij.doseGrid.dimensions(3))-stf(i).isoCenter(3)];

        % isoOffset = [dij.doseGrid.resolution.x dij.doseGrid.resolution.y dij.doseGrid.resolution.z]./2 ...
        %                     +[dij.ctGrid.resolution.x   dij.ctGrid.resolution.y   dij.ctGrid.resolution.z]./2;
        % noClueOffset = [ct.resolution.x ct.resolution.y, ct.resolution.z]./2 + [rem(a(1), a(2)),rem(b(1), b(2)),rem(c(1),c(2))];
        
        isoOffset = [dij.doseGrid.resolution.x dij.doseGrid.resolution.y dij.doseGrid.resolution.z]./2;

        ctDoseGridResolution_x = sort([ct.resolution.x,dij.doseGrid.resolution.x], 'descend');
        ctDoseGridResolution_y = sort([ct.resolution.y,dij.doseGrid.resolution.y], 'descend');
        ctDoseGridResolution_z = sort([ct.resolution.z,dij.doseGrid.resolution.z], 'descend');
        
        noClueOffset = [rem(ctDoseGridResolution_x(1), ctDoseGridResolution_x(2)),...
                        rem(ctDoseGridResolution_y(1), ctDoseGridResolution_y(2)),...
                        rem(ctDoseGridResolution_z(1),ctDoseGridResolution_z(2))];

%        stfFred(i).isoCenter = stfFred(i).isoCenter - isoOffset + noClueOffset;
%        stfFred(i).isoCenter = (stfFred(i).isoCenter - isoOffset + noClueOffset);
        stfFred(i).isoCenter = (-stfFred(i).isoCenter + isoOffset - [ct.resolution.x, ct.resolution.y, ct.resolution.z]);

        nominalEnergies        = unique([stf(i).ray.energy]);
        [~,nominalEnergiesIdx] = intersect([this.machine.data.energy],nominalEnergies);
        
        energyIdxInEmittance   = ismember(emittanceBaseData.energyIndex, nominalEnergiesIdx);

        monteCarloBaseData     = emittanceBaseData.monteCarloData(energyIdxInEmittance);
        
        stfFred(i).nominalEnergies = nominalEnergies;
        stfFred(i).energies        = [monteCarloBaseData.MeanEnergy];

        energySpreadMeV            = [monteCarloBaseData.EnergySpread].*[monteCarloBaseData.MeanEnergy]/100;

        % This energy spread is sigmaE (?) this should match MC2 better
        stfFred(i).energySpread    = energySpreadMeV; %2.355*energySpreadMeV;
        
        stfFred(i).BAMStoIsoDist   = emittanceBaseData.nozzleToIso;%this.machine.meta.BAMStoIsoDist;
    
        [~,eIdx] = intersect(stfFred(i).energies, [this.machine.data.energy]);

        %for the time being use just initial sigma, then need to
        %go through MCemittanceBaseData
        % stfFred(i).FWHMs = [];
        % for j= eIdx'
        %     stfFred(i).FWHMs           = [stfFred(i).FWHMs, [this.machine.data(j).initFocus.sigma(1)]];
        % end
            
        %stfFred(i).FWHMs = 2.355*stfFred(i).FWHMs;
        
        % This is just used for reference and define the field size
        stfFred(i).FWHMs = 2.355*[monteCarloBaseData.SpotSize1x];
        stfFred(i).emittanceX       = [monteCarloBaseData.twissEpsilonX];
        stfFred(i).twissBetaX       = [monteCarloBaseData.twissBetaX];
        stfFred(i).twissAlphaX      = [monteCarloBaseData.twissAlphaX];


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

                if any(stf(i).ray(j).energy == nominalEnergies(k))
                    stfFred(i).energyLayer(k).rayNum   = [stfFred(i).energyLayer(k).rayNum j];
                    stfFred(i).energyLayer(k).bixelNum = [stfFred(i).energyLayer(k).bixelNum ...
                        find(stf(i).ray(j).energy == nominalEnergies(k))];

                    targetX = stf(i).ray(j).targetPoint_bev(1);
                    targetY = stf(i).ray(j).targetPoint_bev(3);

                    sourceX = stf(i).ray(j).rayPos_bev(1);
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

                    stfFred(i).energyLayer(k).rayPosX         = [stfFred(i).energyLayer(k).rayPosX, getPointAtBAMS(targetX,sourceX,distance,stfFred(i).BAMStoIsoDist)];
                    stfFred(i).energyLayer(k).rayPosY         = [stfFred(i).energyLayer(k).rayPosY, getPointAtBAMS(targetY,sourceY,distance,stfFred(i).BAMStoIsoDist)];

                    stfFred(i).energyLayer(k).rayDivX         = [stfFred(i).energyLayer(k).rayDivX, divergenceX];
                    stfFred(i).energyLayer(k).rayDivY         = [stfFred(i).energyLayer(k).rayDivY, divergenceY];
                    
                    if this.calcDoseDirect
                         stfFred(i).energyLayer(k).numOfPrimaries = [stfFred(i).energyLayer(k).numOfPrimaries ...
                                              stf(i).ray(j).weight(stf(i).ray(j).energy == nominalEnergies(k))];
                    else
                         %matRad_cfg.dispWarning('DIJ calculation not yet implemented');
                         stfFred(i).energyLayer(k).numOfPrimaries = [stfFred(i).energyLayer(k).numOfPrimaries ...
                             1];
                    end

                end

            end
        end

        % stfFred(i).nBixels = [];
        % 
        % for k=1:numel(stfFred(i).energies)
        %     stfFred(i).nBixels = [stfFred(i).nBixels, numel(stfFred(i).energyLayer(k).rayNum)];
        % end

        %FRED works in cm
        stfFred(i).isoCenter       = stfFred(i).isoCenter/10;
        stfFred(i).BAMStoIsoDist   = stfFred(i).BAMStoIsoDist/10;
        %stfFred(i).FWHMs           = stfFred(i).FWHMs/10;

        stfFred(i).emittanceX       = stfFred(i).emittanceX./10;
        stfFred(i).twissBetaX       = stfFred(i).twissBetaX./10;

        stfFred(i).totalNumOfBixels = stf(i).totalNumOfBixels;
        for j=1:size(stfFred(i).energies,2)
           stfFred(i).energyLayer(j).rayPosX      = stfFred(i).energyLayer(j).rayPosX/10;
           stfFred(i).energyLayer(j).rayPosY      = stfFred(i).energyLayer(j).rayPosY/10;
           stfFred(i).energyLayer(j).targetPoints = stfFred(i).energyLayer(j).targetPoints/10;
           stfFred(i).energyLayer(j).nBixels      = numel(stfFred(i).energyLayer(j).rayPosX);
           %This is necessary because of the sum(w) at the end of
           %calcDoseForward
           if this.calcDoseDirect
               stfFred(i).energyLayer(j).numOfPrimaries = this.conversionFactor*stfFred(i).energyLayer(j).numOfPrimaries/cumulativeWeights;
           end
        end
    end
    
    
    % %% MC computation and dij filling

    this.writeFredInputAllFiles(stfFred);

    if ~this.exportCalculation
        %should add checks for installation of fred and so on
        matRad_cfg.dispInfo('calling FRED');


        %Need to make this better
        cd(this.MCrunFolder);
        [status,cmdout] = system('fred -f fred.inp','-echo');
        cd(this.FREDrootFolder);
        
        if status==0
            matRad_cfg.dispInfo('done\n');
        end

        if ~this.calcDoseDirect

            %Need to add sanity check for presence of Dij.bin
            %Now working for one field. Need to check what happens with two
            dijFileName = fullfile(this.MCrunFolder, 'out', 'Dij.bin');

            %This introduces a permutation in the collected cube
            dijMatrix = this.readSparseDijBin(dijFileName);

            % Check consistency
            if isequal(size(dijMatrix), [dij.doseGrid.numOfVoxels,dij.totalNumOfBixels])
                %When scoring dij, FRED internaly normalizes to 1
                dij.physicalDose{1} = this.conversionFactor*dijMatrix;
            
            end

        else
            % For the time being, just load the dose and resample
            cube = matRad_readMhd(fullfile(this.MCrunFolder, 'out'), 'Dose.mhd');


            dij.physicalDose{1} = sparse(this.VdoseGrid, ones(numel(this.VdoseGrid),1),cube(this.VdoseGrid), dij.doseGrid.numOfVoxels,1);

            if this.calcLET
                cubeLET = matRad_readMhd(fullfile(this.MCrunFolder, 'out'), 'LETd.mhd');

                dij.mLETDose{1} = sparse(this.VdoseGrid, ones(numel(this.VdoseGrid),1),cubeLET(this.VdoseGrid).*cube(this.VdoseGrid), dij.doseGrid.numOfVoxels,1);
            end

            if this.calcBioDose
                switch this.RBEmodel
                    case 'MCN_RBExD'
                        fName = 'McNamara';

                end

                bioDoseCube = matRad_readMhd(fullfile(this.MCrunFolder, 'out', 'RBE'), ['DoseBio_', fName, '.mhd']);
                dij.BioDose{1} = sparse(this.VdoseGrid, ones(numel(this.VdoseGrid),1),bioDoseCube(this.VdoseGrid), dij.doseGrid.numOfVoxels,1);
            end

        end
    else
        matRad_cfg.dispInfo('All files have been generated');
        dij.physicalDose = [];
    end
    
    this.finalizeDose(ct,cst,stf,dij);
    

    cd(currFolder);

end