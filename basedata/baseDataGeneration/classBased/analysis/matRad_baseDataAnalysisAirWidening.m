classdef matRad_baseDataAnalysisAirWidening < matRad_baseDataAnalysis
    properties
       particleID = 2212; %protons
       data;
       
       
       histogramProperties = struct('binResolution', 0.1, ... %in deg/rad     %this can be turnet to class at a certain point
                                    'binNumber', 100, ...,
                                     'Xcoord', [], ...
                                     'XBins', [], ...
                                     'binResolution_theta', 0.001, ...
                                     'binNumber_theta', 100, ...
                                     'XBins_theta', [], ...
                                     'Xcoord_theta', []); %has to be even
       
 
       initFocus = struct('depths', [], ...
                          'sigma', [], ...
                          'sigmaTheta', [], ...
                          'correlation', [], ...
                          'sigma1', [], ...
                          'sigma2', [], ...
                          'weight', [], ...
                          'sigmaTheta1', [], ...
                          'sigmaTheta2', []);
 
       fitParams = struct('fitPositionSingleGaussian', 1, ...
                          'fitAngularSingleGaussian', 1, ...
                          'fitPositionDoubleGaussian', 1, ...
                          'fitAngularDoubleGaussian',1);
       saveVariables = {'initFocus', 'fitParams', 'histogramProperties'};
       saveNamePrefix = 'airWideningAnalysis';

    end

    methods
        function obj = matRad_baseDataAnalysisAirWidening(mainSimulationClass)
            obj@matRad_baseDataAnalysis(mainSimulationClass);
        end

        function loadPhaseSpace(obj, energyIdx)
            %Here load only the phase space data, combine them into a
            %structure data that has 1 element for each phantom

            fprintf(['Loading Data for Energy ', num2str(obj.simulateEnergies(energyIdx)), '\n']);
            fileFolder = [obj.MCparams.runDirectory, filesep, 'Energy', num2str(obj.simulateEnergies(energyIdx)), filesep, 'Results', filesep, 'PhaseSpace', filesep];
            getAngle = @(angle) acos(angle)-pi/2;
            
            %Build the coordinates of 1 dimensional histogram. This should
            %be turned to a class later
            %Position coordinates
            XBins = linspace(0,obj.histogramProperties.binResolution*obj.histogramProperties.binNumber, obj.histogramProperties.binNumber);
            obj.histogramProperties.XBins = XBins - mean(XBins);
            obj.histogramProperties.Xcoord = obj.histogramProperties.XBins(1:end-1) + obj.histogramProperties.binResolution/2;

            %Angle cooridinates
            XBins = linspace(0,obj.histogramProperties.binResolution_theta*obj.histogramProperties.binNumber_theta, obj.histogramProperties.binNumber_theta);
            obj.histogramProperties.XBins_theta = XBins - mean(XBins);
            obj.histogramProperties.Xcoord_theta = obj.histogramProperties.XBins_theta(1:end-1) + obj.histogramProperties.binResolution_theta/2;
            
            %Loop over phantoms and runs
            for phantomIdx = 1:obj.phantoms.nPhantoms
                data_run = cell(8,1);

                for runIdx = 1:obj.MCparams.nRuns-1
                    if runIdx >2
                        fileName = ['Ps_phantom_', num2str(phantomIdx), '_', num2str(runIdx-1), '.phsp'];
                    else
                        fileName = ['Ps_phantom_', num2str(phantomIdx), '.phsp'];

                    end
                    
                    fID = fopen([fileFolder,fileName]);
                    data_raw = textscan(fID, '%f %f %f %f %f %f %f %f');
                    fclose(fID);

                    particle_Idx = find(data_raw{8}==obj.particleID);
                    if ~isempty(particle_Idx)
                        for k=1:7
                            data_run{k} = [data_run{k}; data_raw{k}(particle_Idx)];
                        end
                    end
                end
                obj.data(phantomIdx).depth = obj.phantoms.depths(phantomIdx);
                obj.data(phantomIdx).posX = data_run{1};
                obj.data(phantomIdx).posY = data_run{2};
                obj.data(phantomIdx).thetaX = getAngle(data_run{4});
                obj.data(phantomIdx).thetaY = getAngle(data_run{5});

                
                obj.data(phantomIdx).histoX      = histcounts(obj.data(phantomIdx).posX, obj.histogramProperties.XBins);
                obj.data(phantomIdx).histoThetaX = histcounts(obj.data(phantomIdx).thetaX, obj.histogramProperties.XBins_theta); 
            end
        end

        function analysis(obj)

            fprintf('Performing analysis...');

            singleEnergy_initFocus = struct('depths', [], ...
                          'sigma', [], ...
                          'sigmaTheta', [], ...
                          'correlation', [], ...
                          'sigma1', [], ...
                          'sigma2', [], ...
                          'weight', [], ...
                          'sigmaTheta1', [], ...
                          'sigmaTheta2', []);

            singleEnergy_initFocus.depths = obj.phantoms.depths;
            for phantomIdx =1:obj.phantoms.nPhantoms

                if obj.fitParams.fitPositionSingleGaussian && (sum(obj.data(phantomIdx).histoX)>0)
                    [initFocus_phantom.sigma] = obj.fitSingleGaussian(obj.histogramProperties.Xcoord,obj.data(phantomIdx).histoX, obj.histogramProperties.XBins);
                else
                    initFocus_phantom.sigma = 0;
                end

                if obj.fitParams.fitAngularSingleGaussian && (sum(obj.data(phantomIdx).histoThetaX)>0)
                    [initFocus_phantom.sigmaTheta] = obj.fitSingleGaussian(obj.histogramProperties.Xcoord_theta,obj.data(phantomIdx).histoThetaX, obj.histogramProperties.XBins_theta);
                else
                    initFocus_phantom.sigmaTheta = 0;
                end

                if obj.fitParams.fitPositionDoubleGaussian && ((sum(obj.data(phantomIdx).histoX)>0))
                    [initFocus_phantom.sigma1, initFocus_phantom.sigma2, initFocus_phantom.weight] = obj.fitDoubleGaussian(obj.histogramProperties.Xcoord,obj.data(phantomIdx).histoX, obj.histogramProperties.XBins);
                else
                    initFocus_phantom.sigma1 = 0;
                    initFocus_phantom.sigma2 = 0;
                    initFocus_phantom.weight = 0;
                end

                if obj.fitParams.fitAngularDoubleGaussian && (sum(obj.data(phantomIdx).histoThetaX)>0)
                    [initFocus_phantom.sigmaTheta1, initFocus_phantom.sigmaTheta2] = obj.fitDoubleGaussian(obj.histogramProperties.Xcoord_theta,obj.data(phantomIdx).histoThetaX,obj.histogramProperties.XBins, initFocus_phantom.weight);
                else
                    initFocus_phantom.sigmaTheta1 = 0;
                    initFocus_phantom.sigmaTheta2 = 0;
                end
               
                if sum(obj.data(phantomIdx).posX>0) && sum(obj.data(phantomIdx).thetaX>0)
                    corr = corrcoef(obj.data(phantomIdx).posX, obj.data(phantomIdx).thetaX);
                    initFocus_phantom.correlation = -corr(1,2);
                else
                    initFocus_phantom.correlation = 0;
                end
               
                singleEnergy_initFocus.sigma = [singleEnergy_initFocus.sigma, initFocus_phantom.sigma];
                singleEnergy_initFocus.sigmaTheta = [singleEnergy_initFocus.sigmaTheta, initFocus_phantom.sigmaTheta];
                singleEnergy_initFocus.correlation = [singleEnergy_initFocus.correlation, initFocus_phantom.correlation];

                singleEnergy_initFocus.sigma1 = [singleEnergy_initFocus.sigma1, initFocus_phantom.sigma1];
                singleEnergy_initFocus.sigma2 = [singleEnergy_initFocus.sigma2, initFocus_phantom.sigma2];
                singleEnergy_initFocus.weight = [singleEnergy_initFocus.weight, initFocus_phantom.weight];
                singleEnergy_initFocus.sigmaTheta1 = [singleEnergy_initFocus.sigmaTheta1, initFocus_phantom.sigmaTheta1];
                singleEnergy_initFocus.sigmaTheta2 = [singleEnergy_initFocus.sigmaTheta2, initFocus_phantom.sigmaTheta2];

            end
            
            obj.store_initFocus(singleEnergy_initFocus);
            
            fprintf('done\n');
        end

        function sigma = fitSingleGaussian(obj,x,y,binEdges)

            if ~iscolumn(x)
                x = x';
            end

            if ~iscolumn(y)
                y=y';
            end

            if ~iscolumn(binEdges)
                binEdges=binEdges';
            end
            
            dx = binEdges(2:end) - binEdges(1:end-1);

            gauss1       = @(p,x) (1./(sqrt(2*pi)*p)) * exp(-x.^2./(2*p.^2));
            objFunc     = fittype(gauss1);
            
            start = 0.1;
            lb = 0;
            ub = 100;
            
            [fitObject1] = fit(x,y./sum(y.*dx), objFunc, 'Start', start, 'Lower', lb, 'Upper', ub);
            
            sigma = fitObject1.p;
        end

        function [sigma1,sigma2,w_out] = fitDoubleGaussian(obj,x,y,binEdges,w)
            if ~iscolumn(x)
                x = x';
            end

            if ~iscolumn(y)
                y=y';
            end

            if ~iscolumn(binEdges)
                binEdges=binEdges';
            end
            
            dx = binEdges(2:end) - binEdges(1:end-1);
            
            gauss1       = @(p,x) (1./(sqrt(2*pi)*p)) * exp(-x.^2./(2*p.^2));
            
            if nargin<5
                gauss2 =  @(w, p1, p2, x) ((1-w) * gauss1(p1,x) + w * gauss1(p2,x));
                objFunc     = fittype(gauss2);

                start = [0.01, 0.1, 0.1];
                lb = [0, 0, 0];
                ub = [1,100,100];
            else
                gauss2    = @(p1, p2, x) ((1-w).* gauss1(p1,x) + w.* gauss1(p2,x));
                objFunc     = fittype(gauss2,'coeff',{'p1','p2'}, 'independent', {'x'});

                start = [0.01, 0.1];
                lb = [0, 0];
                ub = [100,100];
            end
            
            [fitObject] = fit(x,y./sum(y.*dx), objFunc, 'Start', start, 'Lower', lb, 'Upper', ub);
            
            sigma1 = fitObject.p1;
            sigma2 = fitObject.p2;
            if nargin<5
                w_out = fitObject.w;
            else
                w_out = w;
            end
        end

        function store_initFocus(obj, singleEnergy_initFocus)
            if size(obj.initFocus,1) == 1 && ~any(~(structfun(@isempty, obj.initFocus))) 
                 obj.initFocus = singleEnergy_initFocus;
            else
                 obj.initFocus = [obj.initFocus; singleEnergy_initFocus];
            end
        end
    end
end