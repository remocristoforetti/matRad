classdef matRad_dose_analysis < matRad_baseDataAnalysis & matRad_baseDataGeneration_dose
    properties
        fitParams = struct('scoringArea', []);
        fitDoseOutput = struct('depths', [], ...
                          'sigma', [], ...
                          'sigma1', [], ...
                          'sigma2', [], ...
                          'weight', []);
    end


    methods
        function obj = matRad_dose_analysis()
            obj@matRad_baseDataAnalysis();
            obj@matRad_baseDataGeneration_dose();

            obj.saveVariables = {'fitParams', 'fitDoseOutput'};
            
            obj.saveNamePrefix = 'doseAnalysis';
        end

        function data = loadData(obj,energyIdx,~)

            matRad_cfg = MatRad_Config.instance();

            matRad_cfg.dispInfo(['Loading Data for Energy ', num2str(obj.simulateEnergies(energyIdx)), '\n']);

            fileFolder = [obj.MCparams.runDirectory, filesep, 'Energy', num2str(obj.simulateEnergies(energyIdx)), filesep, 'Results', filesep, 'DoseToMedium', filesep];
            data = [];
            for phantomIdx = 1:obj.phantoms.nPhantoms
                data_run = [];

                for runIdx = 1:obj.MCparams.nRuns-1
                    if runIdx >2
                        fileName = ['Dose_phantom_', num2str(phantomIdx), '_', num2str(runIdx-1), '.bin'];
                    else
                        fileName = ['Dose_phantom_', num2str(phantomIdx), '.bin'];

                    end
                    
                    fID = fopen([fileFolder,fileName]);
                    data_raw = fread(fID,obj.phantoms.Rbins(phantomIdx)*obj.phantoms.Zbins(energyIdx,phantomIdx), 'double');
                    data_raw = reshape(data_raw, [obj.phantoms.Rbins(phantomIdx),obj.phantoms.Zbins(energyIdx,phantomIdx)]);
                    fclose(fID);

                    data_run(:,:,runIdx) = data_raw;
                end

                data.(['Phantom', num2str(phantomIdx)]) = mean(data_run,3);
            end
        end

        function analysis(obj,energyIdx,data)
            matRad_cfg = MatRad_Config.instance();
            
            combinedData = obj.combineScorers(data);

            %Check for compatibility of radii
            if length(unique(obj.phantoms.rMax))>1
                matRad_cfg.dispError('Incompatible phantom sizes');
            end
            
            %Define radial coordinates and area
            resolutionR = obj.phantoms.rMax/obj.phantoms.Rbins;
            R = linspace(0,obj.phantoms.rMax(1),obj.phantoms.Rbins(1)+1);
            obj.fitParams(energyIdx).scoringArea = pi*(R(2:end).^2 - R(1:end-1).^2);
            obj.fitParams(energyIdx).r           = R(1:end-1) +  resolutionR/2;

            %Define Zcoordinates
            depths = [];
            shift = 0;
            for phantomIdx = 1:obj.phantoms.nPhantoms
                resolutionZ(phantomIdx) = 2*obj.phantoms.HL(energyIdx,phantomIdx)/obj.phantoms.Zbins(energyIdx,phantomIdx);
                curr_depths = [resolutionZ(phantomIdx)/2:resolutionZ(phantomIdx):2*obj.phantoms.HL(energyIdx,phantomIdx)];

                depths = [depths, curr_depths+shift];
                shift = shift + 2*obj.phantoms.HL(energyIdx,phantomIdx);
            end

            obj.fitParams(energyIdx).depths = depths;

            [PDD, r80, peakPos] = obj.getPDD(combinedData,obj.fitParams(energyIdx).depths, obj.fitParams(energyIdx).scoringArea);



            latProfiles = obj.fitLateralProfile(combinedData,obj.fitParams(energyIdx).r,R, obj.fitParams(energyIdx).scoringArea, obj.energyParams.initFocus.initSigma(energyIdx));

            %Interpolate
            PDD     = interp1(depths, PDD, [0,depths], 'linear', 'extrap')';
            
            %sigma Multi
            weight = interp1(obj.fitParams(energyIdx).depths, latProfiles.w, [0,obj.fitParams(energyIdx).depths], 'linear', 'extrap');
            sigma1 = interp1(obj.fitParams(energyIdx).depths, latProfiles.sigma1, [0, obj.fitParams(energyIdx).depths], 'linear', 'extrap');
            sigma2 = interp1(obj.fitParams(energyIdx).depths, latProfiles.sigma2, [0, obj.fitParams(energyIdx).depths], 'linear', 'extrap');

            probIdx = sigma1 < sigma1(1);
            sigma1(probIdx) = sigma1(1);
            sigma1(sigma1>100*sigma1(1)) = max(sigma1(sigma1<100*sigma1(1)));

            sigma1  = sqrt(sigma1.^2 -sigma1(1).^2);
            sigma2  = sqrt(sigma2.^2 -sigma1(1).^2);


            %sigma
            sigma = interp1(obj.fitParams(energyIdx).depths,latProfiles.sigma,[0, obj.fitParams(energyIdx).depths], 'linear', 'extrap');
            sigma(sigma<sigma(1)) = 0;
            sigma(sigma>100*sigma(1)) = max(sigma(sigma<100*sigma(1)));
            
            sigma(sigma>0) = sqrt(sigma(sigma>0).^2 - sigma(1).^2);

            %depths
            depths = [0,obj.fitParams(energyIdx).depths];

            obj.fitDoseOutput(energyIdx).depths = depths;
            obj.fitDoseOutput(energyIdx).PDD    = PDD;
            obj.fitDoseOutput(energyIdx).peakPos = peakPos; 
            obj.fitDoseOutput(energyIdx).range  = r80;
            obj.fitDoseOutput(energyIdx).sigma  = sigma;
            obj.fitDoseOutput(energyIdx).sigma1 = sigma1;
            obj.fitDoseOutput(energyIdx).sigma2 = sigma2;
            obj.fitDoseOutput(energyIdx).weight = weight;
            obj.fitDoseOutput(energyIdx).samples = latProfiles.samples;
            
            visBool = 0;
            
            if visBool
                figure;
                plot(depths,sigma, '.-');
                grid on;
                grid minor;
                ylim([0,12]);

                figure;
                subplot(1,3,1);
                plot(depths,sigma1, '.-');
                grid on;
                grid minor;
                ylim([0,12]);

                subplot(1,3,2);
                plot(depths,sigma2, '.-');
                grid on;
                grid minor;
                ylim([0,100]);

                subplot(1,3,3);
                plot(depths,weight, '.-');
                grid on;
                grid minor;

                ylim([0,1]);
            end
        end

        function combinedData = combineScorers(obj,data)
            matRad_cfg = MatRad_Config.instance();
            
            combinedData = [];
            currSize = size(data.Phantom1,1);
            for phantomIdx = 1:obj.phantoms.nPhantoms

                if size(data.(['Phantom', num2str(phantomIdx)]),1) == currSize
                    combinedData = [combinedData,data.(['Phantom', num2str(phantomIdx)])];
                else
                    matRad_cfg.dispError('Phantoms have incompatible sizes, define how to combine them');
                end
            end
        end

        function [PDD,r80,peakPos] = getPDD(obj,data,depths,scoringArea)
            
            PDD = (scoringArea*data)./sum(scoringArea);
            
            newDepths = linspace(depths(1),depths(end), 10*length(depths));
            newPDD = interp1(depths,PDD,newDepths);
            
            [~,maxPDD_idx] = max(newPDD);
            PDD_70_idx = find(newPDD > 0.7*max(newPDD),1,'last');

            r80 = interp1(newPDD(maxPDD_idx:PDD_70_idx), newDepths(maxPDD_idx:PDD_70_idx), 0.8*max(newPDD));
            peakPos = newDepths(maxPDD_idx);

            % Convert J/kg to MeV/g, then divide by fluence (#protons/surface) surface = pi*R[mm]^2
            cf = ((6.24150974*sum(scoringArea))*10^(7))/obj.MCparams.nPrimaries;
            
            PDD = PDD.*cf;

            PDD(PDD<0) = 0;

        end

        function latProfiles = fitLateralProfile(obj,data,r,R,scoringArea,initSigma)
            nDepths = size(data,2);

            sigma = [];
            sampleProfile = [];
            for depthIdx = 1:nDepths
                %profile = data(:,depthIdx)./(sum(scoringArea.*data(:,depthIdx)'));
                profile = data(:,depthIdx);
  

                profile = profile.*scoringArea';
                profile = profile./sum(profile);
                profile = profile./scoringArea';
                sampleProfile(:,depthIdx) = profile;
                
                iFit1 = find(profile > 0.0001*max(profile), 1, 'last');
    	        
                if isempty(iFit1), iFit1 = length(profile); end
                
                if r(iFit1) < r(end)
                    profileFit1 = profile(1:iFit1);
                    rdepthFit1 = r(1:iFit1);
                    radialFitMax1 = rdepthFit1(end);
                else
                    profileFit1 = profile;
                    rdepthFit1 = r;
                    radialFitMax1 = r(end);
                end
                
                if depthIdx>1
                    start = sigma(depthIdx-1)-1;
                    lb    = sigma(depthIdx-1)-0.1*sigma(depthIdx-1);
                    ub    = initSigma+100;
                else
                    start = initSigma;
                    lb    = initSigma;
                    ub    = initSigma + 100;
                end


                
                sigma(depthIdx) = obj.fitSingleRadialGaussian(rdepthFit1,profileFit1,R(1:iFit1+1),start,lb,ub);

                if depthIdx>1
                    start = [0.001, sigma1(depthIdx-1), sigma2(depthIdx-1)];
                    lb    = [0.0001, sigma1(depthIdx-1)-0.05*sigma1(depthIdx-1), sigma2(depthIdx-1)-0.3*sigma2(depthIdx-1)];
                    ub    = [1, initSigma+100, sigma2(depthIdx-1)+0.4*sigma2(depthIdx-1)];%initFocus.initSigma+100];
                else
                    start = [0.001, initSigma, initSigma+10];
                    lb    = [0.0001, 0, 0];
                    ub    = [1, initSigma+100, initSigma+100];
                end

                [sigma1(depthIdx), sigma2(depthIdx), w(depthIdx)] = obj.fitDoubleRadialGaussian(rdepthFit1,profileFit1,R(1:iFit1+1),start,lb,ub,[],0);
                
                % if sigma1(depthIdx) > sigma2(depthIdx)
                %     tmp = [sigma1(depthIdx), sigma2(depthIdx), w(depthIdx)];
                %     sigma1(depthIdx) = tmp(2);
                %     sigma2(depthIdx) = tmp(1);
                %     w(depthIdx)      = 1 - tmp(3);
                % end
            end

            latProfiles.sigma  = sigma;
            latProfiles.sigma1 = sigma1;
            latProfiles.sigma2 = sigma2;
            latProfiles.w      = w;

            latProfiles.samples = sampleProfile;
        end

    end
end