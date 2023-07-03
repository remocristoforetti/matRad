classdef matRad_baseDataGeneration < handle
%% properties
    properties
        simulateEnergies
        energyParams;
        MCparams;
        scorerParams;
        phantoms;

        saveDir;
        workingDir;
    end
%% methods
    methods
        function obj = matRad_baseDataGeneration()
            
            matRad_cfg = MatRad_Config.instance();

            obj.energyParams = struct('simulateRanges', [], ...
                                  'simulateEnergySpread', [], ...
                                  'initFocus', []);

            obj.workingDir = [matRad_cfg.matRadRoot, filesep, 'baseData', filesep, 'baseDataGeneration'];

            if ~exist(obj.workingDir, 'dir')
                mkdir(obj.workingDir);
            end
        end

        
        function updateEnergyParams(obj)
            %Here updates energy parameters every time simulateEnergis
            %changes
            obj.energyParams.nEnergies = size(obj.simulateEnergies,1);

            %range
            roughEnergyRange = @(energy) 0.022.*(energy.^1.77);
            obj.energyParams.simulateRanges = arrayfun(@(energy) roughEnergyRange(energy), obj.simulateEnergies);
            
            %energySpread. this LUT is generated with MCemittanceBaseData
            %for protons, starting from generic base data
            energySpreadLUT.energies = [31.728980 36.798565 41.378823 45.597814 49.535352 53.245263 56.765931 60.125900 63.347093 66.446795 69.438927 72.334902 75.144217 77.874882 80.533723 83.126619 85.658677 88.134366 90.557625 92.931948 95.260451 97.545927 99.790893 101.997626 104.168192 106.304476 108.408202 110.480952 112.524181 114.539232 116.527348 118.489682 120.427304 122.341212 124.232337 126.101551 127.949669 129.777457 131.585634 133.374878 135.145827 136.899084 138.635220 140.354773 142.058255 143.746152 145.418923 147.077007 148.720821 150.350762 151.967209 153.570524 155.161052 156.739124 158.305055 159.859149 161.401695 162.932972 164.453247 165.962777 167.461807 168.950576 170.429311 171.898232 173.357551 174.807471 176.248189 177.679896 179.102774 180.517001 181.922747 183.320179 184.709456 186.090733 187.464160 188.829883 190.188042 191.538773 192.882210 194.218479 195.547707 196.870013 198.185516 199.494329 200.796564 202.092328 203.381726 204.664860 205.941829 207.212731 208.477658 209.736703 210.989956 212.237502 213.479427 214.715814 215.946742 217.172292 218.392539 219.607559 220.817425 222.022207 223.221977 224.416802 225.606748 226.791881 227.972264 229.147960 230.319030 231.485532 232.647526 233.805068 234.958214 236.107018];
            energySpreadLUT.energySpread =  [2.866453 3.326473 3.755067 4.153248 4.525126 4.874823 5.205823 5.719472 4.929327 4.424452 4.059499 3.777324 3.549428 3.359572 3.197674 3.057047 2.933063 2.822386 2.722548 2.631672 2.548306 2.471301 2.399739 2.332872 2.270084 2.210865 2.154785 2.101476 2.050629 2.001975 1.955276 1.910324 1.866940 1.824964 1.784246 1.744654 1.706071 1.668393 1.631514 1.595348 1.559807 1.524815 1.490293 1.456176 1.422393 1.388880 1.355579 1.322431 1.289371 1.256348 1.223302 1.190171 1.156899 1.123420 1.089668 1.055573 1.021059 0.986037 0.950419 0.914088 0.876925 0.838788 0.799491 0.758825 0.716526 0.672243 0.625519 0.575731 0.521949 0.462754 0.395668 0.315632 0.208148 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.00000];

            obj.energyParams.simulateEnergySpread = interp1(energySpreadLUT.energies,energySpreadLUT.energySpread, obj.simulateEnergies);

            %initFocus
            obj.energyParams.initFocus = struct('initSigma', [], ...        %initial sigma in mm
                                                'initThetaSigma', [], ...   %initial sigma angula distribution in rad
                                                'correlation', [], ...
                                                'initSigma1', [], ...       %initial sigma1 in mm in case of double source
                                                'initSigma2', [], ...       %initial sigma2 in mm
                                                'initW', [], ...            %initial weight
                                                'initThetaSigma1', [], ...  %initial sigma angula distribution in rad
                                                'initThetaSigma2', []);  %initial sigma angula distribution in rad

            %initFocus need to be provided from the external at the moment
            %Could also add some checks on consistency
        end
        
        function generateTreeDirectory(obj)
            %Generate the directory tree for simulation
            %MainDir
                %Energy1
                    %Files Simulation
                    %Results
                        %Scorer1
                        %Scorer2
           fprintf('Creating Tree directory\n');
           %main dir
           if isempty(obj.MCparams.runDirectory)
               obj.MCparams.rundDirectory = [pwd,filesep,'defaultMCSimulationDirectory'];
           end
           
           if ~exist(obj.MCparams.runDirectory, 'dir')
               mkdir(obj.MCparams.runDirectory);
           end

           
           for energyIdx = 1:obj.energyParams.nEnergies
               energyFolder = [obj.MCparams.runDirectory, filesep, 'Energy', num2str(obj.simulateEnergies(energyIdx))];
               if ~exist(energyFolder, 'dir')
                   mkdir(energyFolder);
               end

               resultFolder = [energyFolder,filesep,'Results'];
               if ~exist(resultFolder, 'dir')
                   mkdir(resultFolder);
               end
               
               for scorerIdx = 1:obj.scorerParams.nScorers
                   scorerFolder = [resultFolder, filesep, obj.scorerParams.scorers{scorerIdx}];
                   if ~exist(scorerFolder, 'dir')
                        mkdir(scorerFolder);
                   end
               end
           end
        end

        function writeRunFiles(obj)
            
            templateFile = fileread(fullfile(obj.MCparams.templateDir,'proton_run.txt'));
            for energyIdx=1:obj.energyParams.nEnergies
                for runIdx=1:obj.MCparams.nRuns
                    fID = fopen(fullfile(obj.MCparams.runDirectory,['Energy',num2str(obj.simulateEnergies(energyIdx))],['proton_run_', num2str(runIdx), '.txt']), 'w');
                    fprintf(fID,templateFile);
                    fprintf(fID, '\n');
                    fprintf(fID,'i:Ts/NumberOfThreads = %u\n',obj.MCparams.numberOfThreads);
                    fprintf(fID,'i:Ts/ShowHistoryCountAtInterval = %d\n', max([1,obj.MCparams.nPrimaries/10]));
                    fprintf(fID,'i:Ts/Seed = %d\n', max([runIdx, runIdx+obj.MCparams.previousRuns]));
                    fprintf(fID,'includeFile = scorers.txt');
                    fclose(fID);
                end
            end
        end

        function writeSimulationParameters(obj)
            
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispWarning('Virtual function to be implemented by the subclass');
            %virtual Function to be implemented
        end

        function writeScorers(obj)

            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispWarning('Virtual function to be implemented by the subclass');

            %virtual Function to be implemented
        end

        function writeBasicFile(obj)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispWarning('Virtual function to be implemented by the subclass');

            %virtual Function to be implemented
        end

        function writeBeamPosition(obj,fID)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispWarning('Virtual function to be implemented by the subclass');

            %virtual Function to be implemented

        end

        function writeSource(obj, fID)
            
            % fprintf(fID, '#Beam Position\n');
            % fprintf(fID, 's:Ge/BeamPosition/Parent      = "World"\n');
            % fprintf(fID, 's:Ge/BeamPosition/Type        = "Group"\n');
            % fprintf(fID, 'd:Ge/BeamPosition/TransX      = 0 m\n');
            % fprintf(fID, 'd:Ge/BeamPosition/TransY      = 0 m\n');
            % fprintf(fID, 'd:Ge/BeamPosition/TransZ      = -1 * Ge/BeamPosition/BAMStoIsoDis m\n');
            % fprintf(fID, 'd:Ge/BeamPosition/RotX        = 0 deg\n');
            % fprintf(fID, 'd:Ge/BeamPosition/RotY        = 0 deg\n');
            % fprintf(fID, 'd:Ge/BeamPosition/RotZ        = 0 deg\n');
            obj.writeBeamPosition(fID);

            if ~(obj.MCparams.doubleSource)
                fprintf(fID, '#Beam Source\n');
                fprintf(fID, 's:So/MySource/Type                    = "Emittance"\n');
                fprintf(fID, 's:So/MySource/Component               = "BeamPosition"\n');
                fprintf(fID, 's:So/MySource/BeamParticle            = "%s"\n', obj.MCparams.sourceParticle);
                fprintf(fID, 's:So/MySource/Distribution            = "BiGaussian"\n');
                fprintf(fID, 'd:So/MySource/BeamEnergy              = Sim/MySource/BeamEnergy MeV\n');
                fprintf(fID, 'u:So/MySource/BeamEnergySpread        = Sim/MySource/BeamEnergySpread\n');
                fprintf(fID, 'd:So/MySource/SigmaX                  = Sim/MySource/SigmaX mm\n');
                fprintf(fID, 'u:So/MySource/SigmaXprime             = Sim/MySource/SigmaThetaX\n');
                fprintf(fID, 'u:So/MySource/CorrelationX            = Sim/MySource/Correlation\n');
                fprintf(fID, 'd:So/MySource/SigmaY                  = So/MySource/SigmaX mm\n');
                fprintf(fID, 'u:So/MySource/SigmaYPrime             = So/MySource/SigmaXprime\n');
                fprintf(fID, 'u:So/MySource/CorrelationY            = So/MySource/CorrelationX\n');
                fprintf(fID, 'i:So/MySource/NumberOfHistoriesInRun  = Sim/MySource/NumberOfHistories\n');
            else
                fprintf(fID, '#Beam Source\n');
                fprintf(fID, 's:So/MySource1/Type                    = "Emittance"\n');
                fprintf(fID, 's:So/MySource1/Component               = "BeamPosition"\n');
                fprintf(fID, 's:So/MySource1/BeamParticle            = "%s"\n', obj.MCparams.sourceParticle);
                fprintf(fID, 's:So/MySource1/Distribution            = "BiGaussian"\n');
                fprintf(fID, 'd:So/MySource1/BeamEnergy              = Sim/MySource/BeamEnergy MeV\n');
                fprintf(fID, 'u:So/MySource1/BeamEnergySpread        = Sim/MySource/BeamEnergySpread\n');
                fprintf(fID, 'd:So/MySource1/SigmaX                  = Sim/MySource/SigmaNarrow mm\n');
                fprintf(fID, 'u:So/MySource1/SigmaXprime             = Sim/MySource/SigmaThetaNarrow\n');
                fprintf(fID, 'u:So/MySource1/CorrelationX            = 0\n');
                fprintf(fID, 'd:So/MySource1/SigmaY                  = So/MySource1/SigmaX mm\n');
                fprintf(fID, 'u:So/MySource1/SigmaYPrime             = So/MySource1/SigmaXprime\n');
                fprintf(fID, 'u:So/MySource1/CorrelationY            = So/MySource1/CorrelationX\n');
                fprintf(fID, 'i:So/MySource1/NumberOfHistoriesInRun  = Sim/MySource/NumberOfHistories1\n');

                fprintf(fID, 's:So/MySource2/Type                    = "Emittance"\n');
                fprintf(fID, 's:So/MySource2/Component               = "BeamPosition"\n');
                fprintf(fID, 's:So/MySource2/BeamParticle            = "%s"\n', obj.MCparams.sourceParticle);
                fprintf(fID, 's:So/MySource2/Distribution            = "BiGaussian"\n');
                fprintf(fID, 'd:So/MySource2/BeamEnergy              = Sim/MySource/BeamEnergy MeV\n');
                fprintf(fID, 'u:So/MySource2/BeamEnergySpread        = Sim/MySource/BeamEnergySpread\n');
                fprintf(fID, 'd:So/MySource2/SigmaX                  = Sim/MySource/SigmaBroad mm\n');
                fprintf(fID, 'u:So/MySource2/SigmaXprime             = Sim/MySource/SigmaThetaBroad\n');
                fprintf(fID, 'u:So/MySource2/CorrelationX            = 0\n');
                fprintf(fID, 'd:So/MySource2/SigmaY                  = So/MySource2/SigmaX mm\n');
                fprintf(fID, 'u:So/MySource2/SigmaYPrime             = So/MySource2/SigmaXprime\n');
                fprintf(fID, 'u:So/MySource2/CorrelationY            = So/MySource2/CorrelationX\n');
                fprintf(fID, 'i:So/MySource2/NumberOfHistoriesInRun  = Sim/MySource/NumberOfHistories2\n');
            end
        end
            
        function saveParameters(obj)
            
            matRad_cfg = MatRad_Config.instance();
            
            variableName = ['baseDataGeneration_', date(), obj.MCparams.sourceParticle, '.mat'];

            saveDirectory = [obj.workingDir,filesep, 'baseDataParameters'];
            if ~exist(saveDirectory, 'dir')
                mkdir(saveDirectory);
            end
            
            obj.saveDir = saveDirectory;
            
            classProperties = properties(obj);
            for propIdx =1:size(classProperties,1)
                saveStr.(classProperties{propIdx}) = obj.(classProperties{propIdx});
            end

            
            save([saveDirectory, filesep, variableName],'saveStr');
        end


        function retriveMainClass(obj, fileName)
               matRad_cfg = MatRad_Config.instance;
            try
                load(fileName, 'saveStr');
            catch
                matRad_cfg.dispError('No file: ',fileName,' found');
            end

            if exist('saveStr', 'var')
                obj.simulateEnergies = saveStr.simulateEnergies;
                obj.energyParams.initFocus = saveStr.energyParams.initFocus;
                obj.MCparams = saveStr.MCparams;
                obj.phantoms = saveStr.phantoms;
                obj.scorerParams = saveStr.scorerParams;
            end
        end
        %% Setters
        function set.simulateEnergies(obj,values)
            obj.simulateEnergies = values;            
            obj.updateEnergyParams();
        end
    end

end