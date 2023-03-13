classdef matRad_MCemittanceBaseData
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % matRad_MCemmitanceBaseData This is the superclass for MonteCarlo base
    % data calculation
    %
    %
    %
    %
    %
    % References
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2020 the matRad development team.
    %
    % This file is part of the matRad project. It is subject to the license
    % terms in the LICENSE file found in the top-level directory of this
    % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
    % of the matRad project, including this file, may be copied, modified,
    % propagated, or distributed except according to the terms contained in the
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        machine                             %matRad base data machine struct
        bdl_path = ''                       %stores path to generated file
        nozzleToIso                         %Nozzle to Isocenter Distance
        smx                                 %Scanning magnet X to isocenter Distance
        smy                                 %Scanning magnet y to isocenter Distance
        monteCarloData                      %MC Phase space data struct
        selectedFocus                       %array containing selected focus indices per energy
        defaultRelativeEnergySpread = 0;    %default energy spread
        matRad_cfg                          %matRad config
        rangeShifters                       %Stores range shifters

        %To force the phase space approximation even if we have the data
        forceSpectrumApproximation  = false; 
        forceEmittanceApproximation = false;
        forceFixedMU                = true;
    end
    
    properties (SetAccess = private)
        stfCompressed   %measure whether function has additional info about
        %the stf
        problemSigma    % = 1, when there was a problem calculating sigma
        energyIndex     %Indices of calculated energies
    end
    
    methods
        function obj = matRad_MCemittanceBaseData(machine,stf)
            %matRad_MCsquareBaseData construct an instance of the MCsquare
            %Base data format using a focus index
            
            %stfCompressed states whether monteCarloData are calculated for
            %all energies (false) or only for energies which exist in given
            %stf. If function is called without stf stfCompressed = false.
            if nargin < 2 || isempty(stf)
                obj.stfCompressed = false;
            else
                obj.stfCompressed = true;
                obj = obj.getRangeShiftersFromStf(stf);
            end
            
            obj.matRad_cfg = MatRad_Config.instance();
                        
            obj.machine = machine;
            obj.problemSigma = false;
            obj.selectedFocus = ones(numel(machine.data),1) * NaN;
            
            if isfield(machine.meta,'BAMStoIsoDist')
                obj.nozzleToIso = machine.meta.BAMStoIsoDist;
            else
                obj.matRad_cfg.dispWarning('No information on BAMS to isocenter distance. Using generic value of 500mm');
                obj.nozzleToIso = 500;
            end
            
            if all(isfield(machine.meta,{'SAD_x','SAD_y'}))
                obj.smx = machine.meta.SAD_x;
                obj.smy = machine.meta.SAD_y;
            elseif isfield(machine.meta,'SAD')
                SAD = machine.meta.SAD;
                obj.smx = SAD;
                obj.smy = SAD;
            else
                obj.matRad_cfg.dispError('No SAD found!');
            end
            
            obj.monteCarloData = [];
            
            %select needed energies and according focus indices by using stf
            if obj.stfCompressed
                tmp = [stf(:).ray];
                plannedEnergies     = [tmp.energy];
                focusIndex          = [tmp.focusIx];
                [~, ind]            = unique(plannedEnergies);
                plannedEnergies     = plannedEnergies(ind);
                focusIndex          = focusIndex(ind);
                [~ ,obj.energyIndex, ~] = intersect([machine.data(:).energy],plannedEnergies);
                
                %if no stf was refered all energies are chosen, while setting
                %the focus index for all energies to preliminary 1
            else
                plannedEnergies = [machine.data(:).energy];
                focusIndex = ones(size(plannedEnergies));
                [~ ,obj.energyIndex, ~] = intersect([machine.data(:).energy],plannedEnergies);
            end
            
            obj.selectedFocus(obj.energyIndex) = focusIndex;
            
            count = 1;
            for i = 1:length(obj.energyIndex)
                ixE = obj.energyIndex(i);
                
                %look up whether MonteCarlo data are already present in
                %machine file , if so do not recalculate
                if isfield(machine.data(ixE), 'energySpectrum') && ~obj.forceSpectrumApproximation
                    energySpectrum = machine.data(ixE).energySpectrum;
                    if isfield(energySpectrum,'type') && strcmp(energySpectrum.type,'gaussian')
                        energyData.NominalEnergy    = ones(1,4) * machine.data(ixE).energy(:);
                        energyData.MeanEnergy       = machine.data(ixE).energySpectrum.mean(:);
                        energyData.EnergySpread     = machine.data(ixE).energySpectrum.sigma(:);
                    else
                        energyData = obj.fitPhaseSpaceForEnergy(ixE);
                    end
                else
                    energyData = obj.fitPhaseSpaceForEnergy(ixE);
                end

                if isfield(machine.data(ixE),'MUdata') && ~obj.forceFixedMU
                    energyData.ProtonsMU = machine.data(ixE).MUdata.numParticlesPerMU;
                else
                    energyData.ProtonsMU = 1e6; %Standard w calibration
                end
                
                if isfield(machine.data(ixE).initFocus,'emittance') && ~obj.forceEmittanceApproximation 
                    data = [];
                    emittance = machine.data(ixE).initFocus.emittance;
                    if ~strcmpi(emittance.type,'bigaussian')
                        matRad_cfg.dispError('Can not handle emittance of type ''%S''!',emittance.type);
                    end

                    nGauss = 1;
                    if isfield(emittance,'weight') 
                        nGauss = length(emittance.weight) + 1;
                    end

                    if nGauss > 2
                        matRad_cfg.dispError('Can not process more than two Gaussians in Emittance parameterization!');
                    end
                    
                    opticsData.Weight1          = 1;
                    opticsData.SpotSize1x       = emittance.sigmaX(1);
                    opticsData.Divergence1x     = emittance.divX(1);
                    opticsData.Correlation1x    = emittance.corrX(1);
                    opticsData.SpotSize1y       = emittance.sigmaY(1);
                    opticsData.Divergence1y     = emittance.divY(1);
                    opticsData.Correlation1y    = emittance.corrY(1);
                    
                    if nGauss == 1
                        opticsData.Weight2          = 0;
                        opticsData.SpotSize2x       = 0;
                        opticsData.Divergence2x     = 0;
                        opticsData.Correlation2x    = 0;
                        opticsData.SpotSize2y       = 0;
                        opticsData.Divergence2y     = 0;
                        opticsData.Correlation2y    = 0;
                    else
                        opticsData.Weight1          = 1 - emittance.weight(1);
                        opticsData.Weight2          = emittance.weight(1);
                        opticsData.SpotSize2x       = emittance.sigmaX(2);
                        opticsData.Divergence2x     = emittance.divX(2);
                        opticsData.Correlation2x    = emittance.corrX(2);
                        opticsData.SpotSize2y       = emittance.sigmaY(2);
                        opticsData.Divergence2y     = emittance.divY(2);
                        opticsData.Correlation2y    = emittance.corrY(2);
                    end

                    %opticsData.FWHMatIso = 2.355 * sigmaNull;
                    opticsData.FWHMatIso = machine.data(ixE).initFocus.SisFWHMAtIso;
                                
                    tmp = energyData;
                    f = fieldnames(opticsData);
                    for a = 1:length(f)
                        tmp.(f{a}) = opticsData.(f{a});
                    end
                    
                    data = [data; tmp];
                else
                    data = [];
                    tmp = energyData;
                    for j = 1:size(machine.data(ixE).initFocus.sigma,1)
                        
                        %                         tmp = energyData;
                        opticsData = obj.fitBeamOpticsForEnergy(ixE, j);
                        
                        f = fieldnames(opticsData);
                        for a = 1:length(f)
                            if j == 1
                                tmp.(f{a}) = opticsData.(f{a});
                            else
                                tmp.(f{a}) = [tmp.(f{a}), opticsData.(f{a})];
                            end
                        end
                        
                        data = tmp;
                    end
                    
                end

                                
                obj.monteCarloData = [obj.monteCarloData, data];
                count = count + 1;
            end
            
            %throw out warning if there was a problem in calculating the
            %width of the Bragg peak in obj.fitBeamOpticsForEnergy
            if obj.problemSigma
                obj.matRad_cfg.dispWarning('Calculation of FWHM of bragg peak in base data not possible! Using simple approximation for energy spread');
            end
        end
        
        
        function mcDataEnergy = fitPhaseSpaceForEnergy(obj,energyIx)
            %function to calculate mean energy and energy spread used by
            %mcSquare for given energy
            
            %Considers air distance from nozzle to phantom surface
            %used in the machine data. 0 means fitted to vacuum simulations
            %with surface at isocenter
            if ~isfield(obj.machine.meta, 'fitAirOffset')
                fitAirOffset = 0;
                %               warning('Could not find fitAirOffset. Using default value (no correction / fit in vacuum).');
            else
                fitAirOffset = obj.machine.meta.fitAirOffset;
            end
            
            %This was changed. This correction is not precise.
            %In order for this correction to be precise, nozzleToIso should
            %be changed to nozzleToSkin = (stf(i).ray(j).SSD + BAM) - SAD;
            %In case of a box phantom with surface at iso SSD = SAD and
            %nozzleToIso = BAM.
            %%%% !!!!!!!!!!!!!!!!!! Correction has nozzleToIso only if want to
            %%%% simulate in air and reproduce baseData range (that was in vaccum)
            %dR = 0.0011 * (obj.nozzleToIso - fitAirOffset);
            %dR = 0.0011 * (fitAirOffset);
            
            i = energyIx;
            
            mcDataEnergy.NominalEnergy = ones(1, size(obj.machine.data(1).initFocus.dist,1)) * obj.machine.data(i).energy;
            
            newDepths = linspace(0,obj.machine.data(i).depths(end),numel(obj.machine.data(i).depths) * 100);
            newDose   = interp1(obj.machine.data(i).depths, obj.machine.data(i).Z, newDepths, 'spline');
            
            %find FWHM w50 of bragg peak and range of 80% does fall off
            [maxV, maxI] = max(newDose);
            [~, r80ind] = min(abs(newDose(maxI:end) - 0.8 * maxV));
            r80ind = r80ind - 1;
            r80 = interp1(newDose(maxI + r80ind - 1:maxI + r80ind + 1), ...
                newDepths(maxI + r80ind - 1:maxI + r80ind + 1), 0.8 * maxV);% ...
            % + obj.machine.data(i).offset + dR;
            
            %Correct r80 with dR
            %r80 = r80 + dR + obj.machine.data(i).offset;
            %r80 = r80 + dR;

            %r80 = r80 + obj.machine.data(i).offset;

            [~, d50rInd] = min(abs(newDose(maxI:end) - 0.5 * maxV));
            d50rInd = d50rInd - 1;
            d50_r = interp1(newDose(maxI + d50rInd - 1:maxI + d50rInd + 1), ...
                newDepths(maxI + d50rInd - 1:maxI + d50rInd + 1), 0.5 * maxV);
            
            if (newDose(1) < 0.5 * maxV)
                [~, d50lInd] = min(abs(newDose(1:maxI) - 0.5*maxV));
                d50_l = interp1(newDose(d50lInd - 1:d50lInd + 1), ...
                    newDepths(d50lInd - 1:d50lInd + 1), 0.5 * maxV);
                w50 = d50_r - d50_l;
                %if width left of peak cannot be determined use r80 as width
            else
                % d50_l = newDepths(maxI);
                w50 = r80;
                obj.problemSigma = true;
            end
            
            %calcualte mean energy used my mcSquare with a formula fitted
            %to TOPAS data
            switch obj.machine.meta.radiationMode
                case 'protons'
                    meanEnergy = @(x) 5.762374661332111e-20 * x^9 - 9.645413625310569e-17 * x^8 + 7.073049219034644e-14 * x^7 ...
                        - 2.992344292008054e-11 * x^6 + 8.104111934547256e-09 * x^5 - 1.477860913846939e-06 * x^4 ...
                        + 1.873625800704108e-04 * x^3 - 1.739424343114980e-02 * x^2 + 1.743224692623838e+00 * x ...
                        + 1.827112816899668e+01;
                    mcDataEnergy.MeanEnergy = ones(1, size(obj.machine.data(1).initFocus.dist,1)) * meanEnergy(r80);
                    
                    %calculate energy straggling using formulae deducted from paper
                    %"An analytical approximation of the Bragg curve for therapeutic
                    %proton beams" by T. Bortfeld et al.
                    totalSigmaSq = ((w50) / 6.14)^2;
                    
                    totalSpreadSq = @(x) 2.713311945114106e-20 * x^9 - 4.267890251195303e-17 * x^8 + 2.879118523083018e-14 * x^7 ...
                        - 1.084418008735459e-11 * x^6 + 2.491796224784373e-09 * x^5 - 3.591462823163767e-07 * x^4 ...
                        + 3.232810400304542e-05 * x^3 - 1.584729282376364e-03 * x^2 + 5.228413840446568e-02 * x ...
                        - 6.547482267336220e-01;
                    
                    % use formula deducted from Bragg Kleeman rule to calcuate
                    % energy straggling given the total sigma and the range
                    % straggling
                    energySpread = (totalSigmaSq - totalSpreadSq(r80)) / (0.022^2 * 1.77^2 * meanEnergy(r80)^(2*1.77-2));
                    energySpread(energySpread < 0) = 0;
                    mcDataEnergy.EnergySpread = ones(1, size(obj.machine.data(1).initFocus.dist,1)) * sqrt(energySpread);
                case 'carbon'
                    % Fit to Range-Energy relationship
                    % Data from "Update to ESTAR, PSTAR, and ASTAR Databases" - ICRU Report 90, 2014
                    % Normalized energy before fit (MeV/u)! Only used ranges [10 350] mm for fit
                    % https://www.nist.gov/system/files/documents/2017/04/26/newstar.pdf
                    meanEnergy = @(x) 11.39 * x^0.628 + 11.24;
                    
                    %%% New fit on HITgantry %%%
                    %meanEnergy = @(x) 10.63 + 13.9*x^0.582 + 0.008584*x^1.422;
                    
                    mcDataEnergy.MeanEnergy = meanEnergy(r80);
                    % reading in a potential given energyspread could go here directly. How would you parse the energyspread
                    % into the function? Through a field in the machine?
                    mcDataEnergy.EnergySpread = obj.defaultRelativeEnergySpread;
                case 'helium'
                    % Fit to Range-Energy relationship
                    % Data from "Update to ESTAR, PSTAR, and ASTAR Databases" - ICRU Report 90, 2014
                    % Normalized energy before fit (MeV/u)! Only used ranges [10 350] mm for fit
                    % https://www.nist.gov/system/files/documents/2017/04/26/newstar.pdf
                    meanEnergy = @(x) 7.57* x.^0.5848 + 3.063;
                    mcDataEnergy.MeanEnergy = meanEnergy(r80);
                    mcDataEnergy.EnergySpread = obj.defaultRelativeEnergySpread;
                otherwise
                    error('not implemented')
            end
        end
        
        
        
        
        function mcDataOptics = fitBeamOpticsForEnergy(obj,energyIx, focusIndex)
            %function to calculate beam optics used by mcSquare for given
            %energy
            
            i = energyIx;
            
            %calculate geometric distances and extrapolate spot size at nozzle
            SAD = obj.machine.meta.SAD;
            z     = -(obj.machine.data(i).initFocus.dist(focusIndex,:) - SAD);
            sigma = obj.machine.data(i).initFocus.sigma(focusIndex,:);
            sigmaSq = sigma.^2;
            
            
            % fitting for either matlab or octave_core_file_limit
            if ~obj.matRad_cfg.isOctave
                
                %fit Courant-Synder equation to data using ipopt, formulae
                %given in mcSquare documentation
                sigmaNull = sqrt(interp1(z,sigmaSq,0));
                
                qRes = @(rho, sigmaT) (sigmaSq -  (sigmaNull^2 - 2*sigmaNull*rho*sigmaT.*z + sigmaT^2.*z.^2));
                
                funcs.objective = @(x) sum(qRes(x(1), x(2)).^2);
                funcs.gradient  = @(x) [  2 * sum(qRes(x(1), x(2)) .* (2 * sigmaNull * x(2) * z));
                    2 * sum(qRes(x(1), x(2)) .* (2 * sigmaNull * x(1) * z  - 2 * x(2) * z.^2))];
                
                options.lb = [-0.99, -Inf];
                options.ub = [ 0.99,  Inf];
                
                options.ipopt.hessian_approximation = 'limited-memory';
                options.ipopt.limited_memory_update_type = 'bfgs';
                options.ipopt.print_level = 1;
                
                start = [0.9; 0.1];
                [result, ~] = ipopt (start, funcs, options);
                rho    = result(1);
                sigmaT = result(2);
                
            else
                
                %fit Courant-Synder equation to data using ipopt, formulae
                %given in mcSquare documentation
                sigmaNull = sqrt(interp1(z,sigmaSq,0));
                
                qRes = @(rho, sigmaT) (sigmaSq -  (sigmaNull^2 - 2*sigmaNull*rho*sigmaT.*z + sigmaT^2.*z.^2));
                
                phi{1} = @(x) sum(qRes(x(1), x(2)).^2);
                phi{2} = @(x) [  2 * sum(qRes(x(1), x(2)) .* (2 * sigmaNull * x(2) * z));
                    2 * sum(qRes(x(1), x(2)) .* (2 * sigmaNull * x(1) * z  - 2 * x(2) * z.^2))];
                
                lb = [-0.99, -Inf];
                ub = [ 0.99,  Inf];
                
                start = [0.9; 0.1];
                [result, ~] = sqp (start, phi, [], [], lb, ub);
                rho    = result(1);
                sigmaT = result(2);
                
            end
            
            %calculate divergence, spotsize and correlation at nozzle
            DivergenceAtNozzle  = sigmaT;
            SpotsizeAtNozzle    = sqrt(sigmaNull^2 - 2 * rho * sigmaNull * sigmaT * obj.nozzleToIso + sigmaT^2 * obj.nozzleToIso^2);
            CorrelationAtNozzle = (rho * sigmaNull - sigmaT * obj.nozzleToIso) / SpotsizeAtNozzle;
            
            
            %save calcuated beam optics data in mcData
            mcDataOptics.ProtonsMU     = 1e6;
            
            mcDataOptics.Weight1       = 1;
            mcDataOptics.SpotSize1x    = SpotsizeAtNozzle;
            mcDataOptics.Divergence1x  = DivergenceAtNozzle;
            mcDataOptics.Correlation1x = CorrelationAtNozzle;
            mcDataOptics.SpotSize1y    = SpotsizeAtNozzle;
            mcDataOptics.Divergence1y  = DivergenceAtNozzle;
            mcDataOptics.Correlation1y = CorrelationAtNozzle;
            
            visBool = false;
            if visBool
                figure, plot(z,sigmaSq,'x');
                zNew = linspace(z(1),z(end),100);
                y = sigmaNull^2 - 2*rho*sigmaNull*sigmaT * zNew + sigmaT^2 * zNew.^2;
                hold on; plot(zNew,y);
            end
            
            
            mcDataOptics.Weight2       = 0;
            mcDataOptics.SpotSize2x    = 0;
            mcDataOptics.Divergence2x  = 0;
            mcDataOptics.Correlation2x = 0;
            mcDataOptics.SpotSize2y    = 0;
            mcDataOptics.Divergence2y  = 0;
            mcDataOptics.Correlation2y = 0;
            mcDataOptics.FWHMatIso = 2.355 * sigmaNull;
        end
        
        
        
        function obj = saveMatradMachine(obj,name)
            %save previously calculated monteCarloData in new baseData file
            %with given name
            
            %             [~ ,energyIndex, ~] = intersect([obj.machine.data(:).energy], [obj.monteCarloData(:).NominalEnergy]);
            
            machineName = [obj.machine.meta.radiationMode, '_', name];
            
            count = 1;
            for i = 1:length(obj.energyIndex)
                
                ixE = obj.energyIndex(i);

                obj.machine.data(ixE).initFocus.emittance.type  = 'bigaussian';
                obj.machine.data(ixE).initFocus.emittance.sigmaX = [obj.monteCarloData(:,count).SpotSize1x obj.monteCarloData(:,count).SpotSize2x];
                obj.machine.data(ixE).initFocus.emittance.sigmaY = [obj.monteCarloData(:,count).SpotSize1y obj.monteCarloData(:,count).SpotSize2y];
                obj.machine.data(ixE).initFocus.emittance.divX = [obj.monteCarloData(:,count).Divergence1x obj.monteCarloData(:,count).Divergence2x];
                obj.machine.data(ixE).initFocus.emittance.divY = [obj.monteCarloData(:,count).Divergence1y obj.monteCarloData(:,count).Divergence2y];
                obj.machine.data(ixE).initFocus.emittance.corrX = [obj.monteCarloData(:,count).Correlation1x obj.monteCarloData(:,count).Correlation2x];
                obj.machine.data(ixE).initFocus.emittance.corrY = [obj.monteCarloData(:,count).Correlation1y obj.monteCarloData(:,count).Correlation2y];
                obj.machine.data(ixE).initFocus.emittance.weight = [obj.monteCarloData(:,count).Weight2]; %Weight one will not be stored explicitly due to normalization
                
                obj.machine.data(ixE).energySpectrum.type  = 'gaussian';
                obj.machine.data(ixE).energySpectrum.mean   = [obj.monteCarloData(:,count).MeanEnergy];
                obj.machine.data(ixE).energySpectrum.sigma = [obj.monteCarloData(:,count).EnergySpread];
                
                
                count = count + 1;
            end
            machine = obj.machine;
            machineFilePath = fullfile(obj.matRad_cfg.matRadRoot,'basedata',[machineName '.mat']);

            save('-v7',machineFilePath,'machine');
            obj.matRad_cfg.dispInfo('Saved Emittance to matRad base data in %s\n',machineFilePath);
        end
    end
    
    methods (Access = protected)
        function obj = getRangeShiftersFromStf(obj,stf)
            allRays = [stf.ray];
            raShis = [allRays.rangeShifter];
            
            [~,ix] =  unique(cell2mat(squeeze(struct2cell(raShis))'),'rows');
            
            raShis = raShis(ix);
            
            ix = [raShis.ID] == 0;
            
            obj.rangeShifters = raShis(~ix);
        end
    end
end

