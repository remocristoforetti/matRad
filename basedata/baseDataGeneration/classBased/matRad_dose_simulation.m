classdef matRad_dose_simulation < matRad_baseDataGeneration_dose
    properties

    end

    
    methods
        function obj = matRad_dose_simulation()
            obj@matRad_baseDataGeneration_dose();
        end

        function writeSimulationParameters(obj, templateFile,energyIdx)

            fID = fopen(fullfile(obj.MCparams.runDirectory,['Energy',num2str(obj.simulateEnergies(energyIdx))],'simulationParameters.txt'), 'w');
            fprintf(fID,templateFile);
            fprintf(fID, '\n');


            for phantomIdx = 1:obj.phantoms.nPhantoms
                phantomName = ['Phantom', num2str(phantomIdx)];
                fprintf(fID, 'd:Sim/%s/HL = %3.3f mm\n', phantomName, obj.phantoms.HL(energyIdx,phantomIdx));
                fprintf(fID, 'u:Sim/%s/ZBins = %3u\n', phantomName, obj.phantoms.Zbins(energyIdx,phantomIdx));
            end
            fprintf(fID, 'd:Sim/MySource/BeamEnergy = %3.3f MeV\n',obj.simulateEnergies(energyIdx));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            energySpread = obj.energyParams.simulateEnergySpread(energyIdx);
            %energySpread = (obj.energyParams.simulateEnergySpread(energyIdx)*100)/obj.simulateEnergies(energyIdx);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf(fID, 'u:Sim/MySource/BeamEnergySpread = %3.3f\n', energySpread);
            if ~obj.MCparams.doubleSource
                fprintf(fID, 'd:Sim/MySource/SigmaX = %3.3f mm\n', obj.energyParams.initFocus.initSigma(energyIdx,1));
                fprintf(fID, 'u:Sim/MySource/SigmaThetaX = %3.4f\n',obj.energyParams.initFocus.initThetaSigma(energyIdx,1));
                fprintf(fID, 'u:Sim/MySource/Correlation = %3.4f\n',obj.energyParams.initFocus.correlation(energyIdx,1));
                fprintf(fID, 'i:Sim/MySource/NumberOfHistories = %u \n',obj.MCparams.nPrimaries);

            else
                nPrimaries1 = floor((1-obj.energyParams.initFocus.initW(energyIdx))*obj.MCparams.nPrimaries);
                nPrimaries2 = floor(obj.energyParams.initFocus.initW(energyIdx)*obj.MCparams.nPrimaries);
                fprintf(fID, 'd:Sim/MySource/SigmaNarrow = %3.3f mm\n', obj.energyParams.initFocus.initSigma1(energyIdx,1));
                fprintf(fID, 'u:Sim/MySource/SigmaThetaNarrow = %3.4f\n',obj.energyParams.initFocus.initThetaSigma1(energyIdx,1));
                fprintf(fID, 'd:Sim/MySource/SigmaBroad = %3.3f mm\n', obj.energyParams.initFocus.initSigma2(energyIdx,1));
                fprintf(fID, 'u:Sim/MySource/SigmaThetaBroad = %3.4f\n',obj.energyParams.initFocus.initThetaSigma2(energyIdx,1));
                fprintf(fID, 'i:Sim/MySource/NumberOfHistories1 = %u \n',nPrimaries1);
                fprintf(fID, 'i:Sim/MySource/NumberOfHistories2 = %u \n',nPrimaries2);
            end

            if any(obj.scorerParams.energyBinned)
                fprintf(fID, 'd:Sim/EMin_protons  = %3.3f MeV\n',obj.scorerParams.Ebinning.EMin);
                fprintf(fID, 'd:Sim/EMax_protons  = %3.3f MeV\n',obj.scorerParams.Ebinning.EMax);
                fprintf(fID, 'i:Sim/EBins_protons  = %u\n',obj.scorerParams.Ebinning.nEBins);           
            end

            fprintf(fID, 'd:Ge/BeamPosition/BAMStoIsoDis = %3.3f mm\n',obj.MCparams.BAMtoISO);

            fclose(fID);
        end


        function writeBasicFile(obj, templateFile)
            fID = fopen(fullfile(obj.MCparams.runDirectory,'proton_basic_dose.txt'), 'w');
            fprintf(fID,templateFile);
            fprintf(fID, '\n');
            for phantomIdx=1:obj.phantoms.nPhantoms
               phantomName = ['Phantom', num2str(phantomIdx)];

               fprintf(fID, '\n');
               fprintf(fID, 's:Ge/%s/Parent            = "MyPhantom"\n', phantomName);
               fprintf(fID, 's:Ge/%s/Type              = "TsCylinder"\n', phantomName);
               fprintf(fID, 's:Ge/%s/Material          = "%s"\n', phantomName, obj.phantoms.material);
               fprintf(fID, 'd:Ge/%s/RMin              = 0.0 cm\n', phantomName);
               fprintf(fID, 'd:Ge/%s/RMax              = %3.3f mm\n', phantomName, obj.phantoms.rMax(phantomIdx));
               fprintf(fID, 'd:Ge/%s/HL                = Sim/%s/HL mm\n', phantomName, phantomName);
               fprintf(fID, 'd:Ge/%s/SPhi              = 0.0 deg\n', phantomName);
               fprintf(fID, 'd:Ge/%s/DPhi              = 360.0 deg\n', phantomName);
               fprintf(fID, 'i:Ge/%s/ZBins             = Sim/%s/ZBins\n', phantomName,phantomName);
               fprintf(fID, 'i:Ge/%s/RBins             = %u\n', phantomName, obj.phantoms.Rbins(phantomIdx));
               
               if phantomIdx == 1
                   fprintf(fID, 'd:Ge/%s/TransZ            = Ge/%s/HL mm\n', phantomName, phantomName);
                   fprintf(fID, 'd:Ge/%s/Edge              = 2 * Ge/%s/HL mm\n', phantomName,phantomName);
               else
                   fprintf(fID, 'd:Ge/%s/TransZ            = Ge/%s/Edge + Ge/%s/HL mm\n', phantomName,['Phantom', num2str(phantomIdx-1)], phantomName);
                   fprintf(fID, 'd:Ge/%s/Edge              = Ge/%s/TransZ + Ge/%s/HL mm\n', phantomName,phantomName,phantomName);
               end


               fprintf(fID, '\n');
               fprintf(fID, '#################################################################################\n');

            end
           obj.writeSource(fID);
           fclose(fID);
        end

        function writeBeamPosition(obj,fID)

            fprintf(fID, '#Beam Position\n');
            fprintf(fID, 's:Ge/BeamPosition/Parent      = "World"\n');
            fprintf(fID, 's:Ge/BeamPosition/Type        = "Group"\n');
            fprintf(fID, 'd:Ge/BeamPosition/TransX      = 0 m\n');
            fprintf(fID, 'd:Ge/BeamPosition/TransY      = 0 m\n');
            fprintf(fID, 'd:Ge/BeamPosition/TransZ      = %d mm\n', obj.phantoms.sourcePosition);
            fprintf(fID, 'd:Ge/BeamPosition/RotX        = 0 deg\n');
            fprintf(fID, 'd:Ge/BeamPosition/RotY        = 0 deg\n');
            fprintf(fID, 'd:Ge/BeamPosition/RotZ        = 0 deg\n');

        end

        function writeScorerIncluder (obj, templateFile, energyIdx)
            fID = fopen(fullfile(obj.MCparams.runDirectory,['Energy',num2str(obj.simulateEnergies(energyIdx))],'scorers.txt'), 'w');
            fprintf(fID,templateFile);
            fprintf(fID, '\n');

            for scorerIdx = 1:obj.scorerParams.nScorers
                scorerName = obj.scorers{scorerIdx};
                %For the time being this is one, then will become ion-dependend

                if  obj.scorerParams.filteredScorer(scorerIdx)
                    scorerName = [scorerName, '_', obj.scorerParams.ions{1}];
                end
           
                fprintf(fID, 'includeFile = ../scorer_%s.txt\n', scorerName);
            end
            fclose(fID);
        end

        function writeScorers(obj,templateFile,scorerIdx)
            scorerName = obj.scorers{scorerIdx};
            %For the time being this is one, then will become ion-dependend
            if  obj.scorerParams.filteredScorer(scorerIdx)
                scorerName = [scorerName, '_', obj.scorerParams.ions{1}];
            end
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
                fprintf(fID, 's:Sc/phantom_%u/%s/OutputFile                   = "./Results/%s/%s_phantom_%u"\n',phantomIdx,scorerName,scorerName,scorerName,phantomIdx);
                if obj.scorerParams.filteredScorer(scorerIdx)
                    fprintf(fID, 'i:Sc/phantom_%u/%s/OnlyIncludeParticlesOfAtomicNumber                   = 1\n',phantomIdx,scorerName);
                end

                if strcmp(obj.scorers{scorerIdx}, 'ProtonLET')
                    fprintf(fID, 'd:Sc/phantom_%u/%s/MaxScoredLET = 100 MeV/mm/(g/cm3) # default 100 MeV/mm/(g/cm3)\n',phantomIdx,scorerName);
                end
                if obj.scorerParams.energyBinned(scorerIdx)
                    fprintf(fID, 'd:Sc/phantom_%u/%s/EBinMax                   = Sim/EMax_protons MeV\n',phantomIdx,scorerName);
                    fprintf(fID, 'd:Sc/phantom_%u/%s/EBinMin                   = Sim/EMin_protons MeV\n',phantomIdx,scorerName);
                    fprintf(fID, 'i:Sc/phantom_%u/%s/EBins                     = Sim/EBins_protons\n',phantomIdx,scorerName);
                    fprintf(fID, 's:Sc/phantom_%u/%s/EBinEnergy                = "PreStep"\n',phantomIdx,scorerName);
                end
                fprintf(fID, '\n');
            end
            fclose(fID);
        end

        function interpInitFocus(obj,initFocusStruct)
            for energyIdx = 1:obj.energyParams.nEnergies
                obj.energyParams.initFocus.initSigma(energyIdx,1)         = interp1(initFocusStruct(energyIdx).depths, initFocusStruct(energyIdx).sigma, obj.MCparams.BAMtoISO);
                obj.energyParams.initFocus.initThetaSigma(energyIdx,1)    = interp1(initFocusStruct(energyIdx).depths, initFocusStruct(energyIdx).sigmaTheta, obj.MCparams.BAMtoISO);
                obj.energyParams.initFocus.correlation(energyIdx,1)       = interp1(initFocusStruct(energyIdx).depths, initFocusStruct(energyIdx).correlation, obj.MCparams.BAMtoISO);
                obj.energyParams.initFocus.initSigma1(energyIdx,1)        = interp1(initFocusStruct(energyIdx).depths, initFocusStruct(energyIdx).sigma1, obj.MCparams.BAMtoISO);
                obj.energyParams.initFocus.initSigma2(energyIdx,1)        = interp1(initFocusStruct(energyIdx).depths, initFocusStruct(energyIdx).sigma2, obj.MCparams.BAMtoISO);
                obj.energyParams.initFocus.initW(energyIdx,1)             = interp1(initFocusStruct(energyIdx).depths, initFocusStruct(energyIdx).weight, obj.MCparams.BAMtoISO);
                obj.energyParams.initFocus.initThetaSigma1(energyIdx,1)   = interp1(initFocusStruct(energyIdx).depths, initFocusStruct(energyIdx).sigmaTheta1, obj.MCparams.BAMtoISO);
                obj.energyParams.initFocus.initThetaSigma2(energyIdx,1)   = interp1(initFocusStruct(energyIdx).depths, initFocusStruct(energyIdx).sigmaTheta2, obj.MCparams.BAMtoISO);
           end
        end

    end
end