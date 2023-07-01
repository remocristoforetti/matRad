classdef matRad_airWidening_simulation < matRad_baseDataGeneration_airWidening
    properties

    end


    methods
        function obj = matRad_airWidening_simulation
            obj@matRad_baseDataGeneration_airWidening();
        end

        function writeSimulationParameters(obj)
            templateFile = fileread(fullfile(obj.MCparams.templateDir,'simulationParameters.txt'));
            for energyIdx=1:obj.energyParams.nEnergies
                    fID = fopen(fullfile(obj.MCparams.runDirectory,['Energy',num2str(obj.simulateEnergies(energyIdx))],'simulationParameters.txt'), 'w');
                    fprintf(fID,templateFile);
                    fprintf(fID, '\n');

                    fprintf(fID, 'd:Sim/HL = %3.3f mm\n', obj.phantoms.HL);
                    fprintf(fID, 'd:Sim/MySource/BeamEnergy = %3.3f MeV\n',obj.simulateEnergies(energyIdx));
                    fprintf(fID, 'u:Sim/MySource/BeamEnergySpread = %3.3f\n', obj.energyParams.simulateEnergySpread(energyIdx));
                    if ~obj.MCparams.doubleSource
                        fprintf(fID, 'd:Sim/MySource/SigmaX = %3.3f mm\n', obj.energyParams.initFocus.initSigma(energyIdx));
                        fprintf(fID, 'u:Sim/MySource/SigmaThetaX = %3.4f\n',obj.energyParams.initFocus.initThetaSigma(energyIdx));
                        fprintf(fID, 'u:Sim/MySource/Correlation = %3.4f\n',obj.energyParams.initFocus.correlation(energyIdx));
                        fprintf(fID, 'i:Sim/MySource/NumberOfHistories = %u \n',obj.MCparams.nPrimaries);

                    else
                        nPrimaries1 = floor((1-obj.energyParams.initFocus.initW)*obj.MCparams.nPrimaries);
                        nPrimaries2 = floor(obj.energyParams.initFocus.initW*obj.MCparams.nPrimaries);
                        fprintf(fID, 'd:Sim/MySource/SigmaNarrow = %3.3f mm\n', obj.energyParams.initFocus.initSigma1(energyIdx));
                        fprintf(fID, 'u:Sim/MySource/SigmaThetaNarrow = %3.3f\n',obj.energyParams.initFocus.initThetaSigma1(energyIdx));
                        fprintf(fID, 'd:Sim/MySource/SigmaBroad = %3.3f mm\n', obj.energyParams.initFocus.initSigma2(energyIdx));
                        fprintf(fID, 'u:Sim/MySource/SigmaThetaBroad = %3.3f\n',obj.energyParams.initFocus.initThetaSigma2(energyIdx));
                        fprintf(fID, 'i:Sim/MySource/NumberOfHistories1 = %u \n',nPrimaries1);
                        fprintf(fID, 'i:Sim/MySource/NumberOfHistories2 = %u \n',nPrimaries2);
                    end

                    fprintf(fID, 'd:Ge/BeamPosition/BAMStoIsoDis = %3.3f mm\n',obj.MCparams.BAMtoISO);

                    fclose(fID);
            end
        end

        function writeScorers(obj)
            templateFile = fileread(fullfile(obj.MCparams.templateDir,'scorer.txt'));
            for energyIdx=1:obj.energyParams.nEnergies
                fID = fopen(fullfile(obj.MCparams.runDirectory,['Energy',num2str(obj.simulateEnergies(energyIdx))],'scorers.txt'), 'w');
                fprintf(fID,templateFile);
                fprintf(fID, '\n');
                for scorerIdx = 1:obj.scorerParams.nScorers
                    fprintf(fID, 'includeFile = ../scorer_%s.txt\n', obj.scorerParams.scorers{scorerIdx});
                end
                fclose(fID);
            end

            for scorerIdx = 1:obj.scorerParams.nScorers
                scorerName = obj.scorerParams.scorers{scorerIdx};
                templateFile = fileread(fullfile(obj.MCparams.templateDir,['scorer_',scorerName,'.txt']));
                fID = fopen(fullfile(obj.MCparams.runDirectory,['scorer_',scorerName,'.txt']), 'w');
                fprintf(fID,templateFile);
                fprintf(fID, '\n');
                for phantomIdx = 1:obj.phantoms.nPhantoms
                    fprintf(fID, '\n');
                    fprintf(fID, 's:Sc/phantom_%u/%s/Quantity                     = Sim/ScoredQuantity_%s\n',phantomIdx,scorerName,scorerName);
                    fprintf(fID, 's:Sc/phantom_%u/%s/Component                    = "Phantom%u"\n',phantomIdx,scorerName,phantomIdx);
                    fprintf(fID, 'b:Sc/phantom_%u/%s/OutputToConsole              = "False"\n',phantomIdx,scorerName);
                    fprintf(fID, 's:Sc/phantom_%u/%s/IfOutputFileAlreadyExists    = "Increment"\n',phantomIdx,scorerName);
                    fprintf(fID, 's:Sc/phantom_%u/%s/OutputType                   = Sim/OutputType_%s\n',phantomIdx,scorerName, scorerName);
                    fprintf(fID, 's:Sc/phantom_%u/%s/OutputFile                   = "./Results/%s/Ps_phantom_%u"\n',phantomIdx,scorerName,scorerName,phantomIdx);
                    fprintf(fID, 's:Sc/phantom_%u/%s/Surface                      = "Phantom%u/ZMinusSurface"',phantomIdx,scorerName,phantomIdx);
                    fprintf(fID, '\n');
                end
                fclose(fID);
            end
        end
        
        function writeBasicFile(obj)
            templateFile = fileread(fullfile(obj.MCparams.templateDir,'proton_basic_air_widening.txt'));
            fID = fopen(fullfile(obj.MCparams.runDirectory,'proton_basic_air_widening.txt'), 'w');
            fprintf(fID,templateFile);
            fprintf(fID, '\n');
            for phantomIdx=1:obj.phantoms.nPhantoms
               phantomName = ['Phantom', num2str(phantomIdx)];

               fprintf(fID, '\n');
               fprintf(fID, 's:Ge/%s/Parent            = "MyPhantom"\n', phantomName);
               fprintf(fID, 's:Ge/%s/Type              = "TsCylinder"\n', phantomName);
               fprintf(fID, 's:Ge/%s/Material          = "Air"\n', phantomName);
               fprintf(fID, 'd:Ge/%s/RMin              = 0.0 cm\n', phantomName);
               fprintf(fID, 'd:Ge/%s/RMax              = %3.3f mm\n', phantomName, obj.phantoms.rMax);
               fprintf(fID, 'd:Ge/%s/HL                = %f mm\n', phantomName, obj.phantoms.HL);
               fprintf(fID, 'd:Ge/%s/SPhi              = 0.0 deg\n', phantomName);
               fprintf(fID, 'd:Ge/%s/DPhi              = 360.0 deg\n', phantomName);
               fprintf(fID, 'i:Ge/%s/ZBins             = 1\n', phantomName);
               fprintf(fID, 'i:Ge/%s/RBins             = 1\n', phantomName);
               fprintf(fID, 'd:Ge/%s/TransZ            = %3.3f mm\n', phantomName, -(obj.MCparams.BAMtoISO - obj.phantoms.depths(phantomIdx) - 2*obj.phantoms.HL));
               fprintf(fID, '\n');
               fprintf(fID, '#################################################################################\n');

            end
           obj.writeSource(fID);
           fclose(fID);
        end


        % function retriveMainClass(obj,fileName)
        % 
        %     matRad_cfg = MatRad_Config.instance;
        %     try
        %         load(fileName, 'saveStr');
        %     catch
        %         matRad_cfg.dispError('No file: ',fileName,' found');
        %     end
        % 
        %     if exist('saveStr', 'var')
        %         obj.simulateEnergies = saveStr.simulateEnergies;
        %         obj.energyParams.initFocus = saveStr.energyParams.initFocus;
        %         obj.MCparams = saveStr.MCparams;
        %         obj.phantoms = saveStr.phantoms;
        %         obj.scorerParams = saveStr.scorerParams;
        %     end
        % end
    end
end