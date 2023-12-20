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
    fred_cfg = MatRad_FREDConfig.instance();
    %matRad_cfg.dispError('Still need to implement this function');
    
    currFolder = pwd;

    % cd to FRED folder (necessary ?)

    cd(fred_cfg.FREDrootFolder);

    %Now we can run initDoseCalc as usual
    dij = this.initDoseCalc(ct,cst,stf);

     for s = 1:dij.numOfScenarios
         HUcube{s} =  matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y',  dij.ctGrid.z,ct.cubeHU{s}, ...
                                     dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'linear');
        %HUcube{s} =  ct.cubeHU{s};
     end

    %check for presence of HUcorrection file
    if ~exist(fullfile(fred_cfg.FREDrootFolder, 'HUmaterialConversionTables'),'dir')
        this.useInternalHUConversion = true;
        if any(HUcube{1}(:)>fred_cfg.hLutLimits(2)) || any(HUcube{1}(:)<fred_cfg.hLutLimits(1))
            matRad_cfg.dispWarning('HU outside of boundaries');
            this.HUclamping = true;
        end

    end
   

    %Write the directory tree necessary for the simulation
    this.writeTreeDirectory();
    
    %patientFileName = fullfile(fred_cfg.regionsFolder, fred_cfg.patientFilename);
    cd(fred_cfg.regionsFolder);
    matRad_writeMhd(HUcube{1},[this.doseGrid.resolution.x, this.doseGrid.resolution.y, this.doseGrid.resolution.z], fred_cfg.patientFilename, 'MET_SHORT');
    cd(fred_cfg.FREDrootFolder);

    getPointAtBAMS = @(target,source,distance,BAMStoIso) (target -source)*(-BAMStoIso)/distance + source;%(target  + source*(BAMStoIso - distance))/distance;
          
    % Loop over the stf to rearrange data
    counter = 0;

    for i = 1:length(stf)

        stfFred(i).gantryAngle     = stf(i).gantryAngle;
        stfFred(i).couchAngle      = stf(i).couchAngle;

        %This is Topas like definition of isocenter
        stfFred(i).isoCenter       = -[0.5*ct.resolution.x*(ct.cubeDim(2)+1)-stf(i).isoCenter(1),...
                                      0.5*ct.resolution.y*(ct.cubeDim(1)+1)-stf(i).isoCenter(2),...
                                      0.5*ct.resolution.z*(ct.cubeDim(3)+1)-stf(i).isoCenter(3)];
        
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
                                              stf(i).ray(j).weight(stf(i).ray(j).energy == stfFred(i).energies(k))];
                    else
                         matRad_cfg.dispWarning('DIJ calculation not yet implemented');
                         stfFred(i).energyLayer(k).numOfPrimaries = [stfFred(i).energyLayer(k).numOfPrimaries ...
                             MCsquareConfig.Num_Primaries];
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
        stfFred(i).FWHMs           = stfFred(i).FWHMs/10;
         
        for j=1:size(stfFred(i).energies,2)
           stfFred(i).energyLayer(j).rayPosX      = stfFred(i).energyLayer(j).rayPosX/10;
           stfFred(i).energyLayer(j).rayPosY      = stfFred(i).energyLayer(j).rayPosY/10;
           stfFred(i).energyLayer(j).targetPoints = stfFred(i).energyLayer(j).targetPoints/10;
           stfFred(i).energyLayer(j).nBixels      = numel(stfFred(i).energyLayer(j).rayPosX);
        end
    end
     
    % %% MC computation and dij filling

    this.writeFredInputAllFiles(stfFred);

    if ~this.exportCalculation
        %should add checks for installation of fred and so on
        matRad_cfg.dispInfo('calling FRED');


        %Need to make this better
        cd(fred_cfg.MCrunFolder);
        [status,cmdout] = system('fred -f fred.inp','-echo');
        cd(fred_cfg.FREDrootFolder);
        
        if status==0
            matRad_cfg.dispInfo('done\n');
        end
    else
        matRad_cfg.dispInfo('All files have been generated');
        dij.physicalDose = [];
    end
    
    this.finalizeDose(ct,cst,stf,dij);
    
    cd(currFolder);
end
