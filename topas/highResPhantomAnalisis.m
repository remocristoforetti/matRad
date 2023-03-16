classdef highResPhantomAnalisis < handle
    properties
        nPhantoms;
        wDir;
        readPDD;
        
        readDSSpectra;
        readSpectra;
        readSPeDSpectra;
        readProtonLET;
        depths;
        phantoms;
        PDD;
        Phi;
        edPhi;

        ProtonLET;

        DSPhi;
        run;
        integrationSurface;
        ions;
        ionsZ;
        downSamplingPhantoms;
        newDepths;

        newPDD;
        newPhi;
        newedPhi;
        newEnergies;

    end

    methods
        function obj = highResPhantomAnalisis()

        end

        function importRawData(obj)
            depths = [];
            for k=1:obj.nPhantoms
                nBins = ceil(obj.phantoms(k).dimension(1)/obj.phantoms(k).resolution(1));
                currDepths = linspace(0,obj.phantoms(k).dimension(1)-obj.phantoms(k).resolution(1),nBins) + obj.phantoms(k).resolution(1)/2 + obj.phantoms(k).positionDepth(k);
                %currDepths = [0:obj.phantoms(k).resolution(1):obj.phantoms(k).dimension(1)-obj.phantoms(k).resolution(1)] + obj.phantoms(k).resolution(1)/2 + obj.phantoms(k).positionDepth(k);
                depths = [depths,currDepths];
            end
            obj.depths = depths;

            PDD = [];

            if obj.readPDD
                scorersDir = [obj.wDir,filesep, 'Output',filesep, 'PDD', filesep, 'Phantom'];
                data = [];
                for k = 1:obj.nPhantoms
                    dataFilename = strcat(scorersDir, num2str(k), filesep, 'physicalDose_Run_', compose('%04d', obj.run), '.dcm');

                    data_dicom = double(squeeze(dicomread(dataFilename{1})));
                    data_dicom_info = dicominfo(dataFilename{1});
                    data_dicom_permute = permute(data_dicom, [2,3,1]);
                    data_dicom = squeeze(sum(data_dicom_permute(obj.integrationSurface, obj.integrationSurface, :), [1,2])).*data_dicom_info.DoseGridScaling;
                    PDD = [PDD; data_dicom];
                end
                obj.PDD = PDD;
            end

            if obj.readDSSpectra
                scorersDir = [obj.wDir,filesep, 'Output',filesep, 'DS', filesep, 'Phantom'];
                data = [];

                for m = 1:size(obj.ionsZ,2)
                    Phi = [];
                    edPhi = [];
                    for k = 1:obj.nPhantoms
                        dataFilename = strcat(scorersDir, num2str(k), filesep, 'Ion_', num2str(obj.ionsZ(m)),filesep, 'S_',obj.ions{m}, '_Run_', compose('%04d', obj.run), '.dcm');
                        data_dicom = double(squeeze(dicomread(dataFilename{1})));
                        data_dicom_info = dicominfo(dataFilename{1});
                        data_dicom_permute = permute(data_dicom, [2,3,1]);
                        data_dicom = data_dicom_permute(2:end,2:end,:).*data_dicom_info.DoseGridScaling./obj.phantoms(k).resolution(1); %This is corrected for phantom depth resolution. In the end, computation of alpha/beta does not depend on normalization
                    
                        Ed = obj.phantoms(k).EParam(m).Ed;
                        Ed = repmat(Ed,1,1,size(data_dicom,3));
                        
                        Phi = [Phi, sparse(squeeze(sum(data_dicom, [1])))];
                        edPhi = [edPhi, sparse(squeeze(pagemtimes(Ed,data_dicom)))];
                    
                    end
                    obj.DSPhi{m} = Phi;

                    obj.edPhi{m} = edPhi;
                end

            end

            %These are the energy spectra, count for number of ions, not
            %number of events
            if obj.readSpectra
                scorersDir = [obj.wDir,filesep, 'Output',filesep, 'SP', filesep, 'Phantom'];
                data = [];

                for m = 1:size(obj.ionsZ,2)
                    Phi = [];

                    for k = 1:obj.nPhantoms
                        dataFilename = strcat(scorersDir, num2str(k), filesep, 'Ion_', num2str(obj.ionsZ(m)),filesep, 'SP_',obj.ions{m}, '_Run_', compose('%04d', obj.run), '.csv');
                        data_csv = csvread(dataFilename{1},11);
                        %This is fluence, (1/mm^2), su multiply by phantom
                        %lateral surface (scoring is performed over whole phantom surface)
                        data_csv = data_csv(:,2:end-2)'.*prod(obj.phantoms(k).dimension([2,3]));
                        
                        Phi = [Phi, sparse(squeeze(data_csv))];
                    end
                    obj.Phi{m} = Phi;

                end

            end

            if obj.readSPeDSpectra

                scorersDir = [obj.wDir,filesep, 'Output',filesep, 'SPeD', filesep, 'Phantom'];
                data = [];

                for m = 1:size(obj.ionsZ,2)
                    Phi = [];

                    for k = 1:obj.nPhantoms
                        dataFilename = strcat(scorersDir, num2str(k), filesep, 'Ion_', num2str(obj.ionsZ(m)),filesep, 'SPeD_',obj.ions{m}, '_Run_', compose('%04d', obj.run), '.dcm');
                        data_dicom = double(squeeze(dicomread(dataFilename{1})));
                        data_dicom_info = dicominfo(dataFilename{1});
                        data_dicom = data_dicom*data_dicom_info.DoseGridScaling./obj.phantoms(k).resolution(1);

                        Phi = [Phi, sparse(squeeze(data_dicom(:,2:end)'))];
                    end
                    obj.edPhi{m} = Phi;

                end

            end

             if obj.readProtonLET
                scorersDir = [obj.wDir,filesep, 'Output',filesep, 'ProtonLET', filesep, 'Phantom'];
                data = [];
                LET = [];

                for k = 1:obj.nPhantoms

                    dataFilename = strcat(scorersDir, num2str(k), filesep, obj.ions{1}(1:end-1), 'LET_Run_', compose('%04d', obj.run), '.dcm');
                    data_dicom = double(squeeze(dicomread(dataFilename{1})));
                    data_dicom_info = dicominfo(dataFilename{1});
                    data_dicom = data_dicom*data_dicom_info.DoseGridScaling;

                    data_dicom_permute = permute(data_dicom, [2,3,1]);
                    data_dicom = squeeze(sum(data_dicom_permute, [1,2]))';
                    LET = [LET, squeeze(data_dicom)];
                end
                    obj.ProtonLET = LET;
             end

             obj.rawDataDownsampling();
        end




        function rawDataDownsampling(obj)

     
            res = obj.downSamplingPhantoms.resolutions;
            nRegions = size(res,2);
            highResWindowWidth = obj.downSamplingPhantoms.highResWindowWidth;

            %Start from PDD
            %get r80
            newDepths = linspace(0,obj.depths(end),numel(obj.depths) * 100);
            newDose   = interp1(obj.depths, obj.PDD, newDepths, 'spline');
            [maxV, maxI] = max(newDose);
            [~, r80ind] = min(abs(newDose(maxI:end) - 0.8 * maxV));
            r80ind = r80ind - 1;
            r80 = interp1(newDose(maxI + r80ind - 1:maxI + r80ind + 1), ...
                              newDepths(maxI + r80ind - 1:maxI + r80ind + 1), 0.8 * maxV);

            %Define high res window
            highResWindow = r80 - [highResWindowWidth,-highResWindowWidth]./2;
            nBins = ceil(highResWindowWidth/res(2));

            newHighResDepths = linspace(highResWindow(1),highResWindow(2),nBins);
            newHighResPDD = interp1(obj.depths, obj.PDD, newHighResDepths, 'spline');

            %Define entrance window
            entranceWindowWidth = [0,highResWindow(1)];
            nBins = ceil(entranceWindowWidth(2)/res(1));

            newEntranceDepths = linspace(0,entranceWindowWidth(2),nBins+1);
            newEntranceDepths = newEntranceDepths(1:end-1);
            newEntrancePDD = interp1([0,obj.depths], [obj.PDD(1); obj.PDD], newEntranceDepths,'spline');

            %Define distal window
            distalWindowWidth = [highResWindow(2),obj.depths(end)];
            nBins = ceil((distalWindowWidth(2)-distalWindowWidth(1))/res(3));

            newDistalDepths = linspace(distalWindowWidth(1),distalWindowWidth(2),nBins+1);
            newDistalDepths = newDistalDepths(2:end);
            newDistalPDD = interp1(obj.depths, obj.PDD, newDistalDepths, 'spline');

            obj.newDepths = [newEntranceDepths,newHighResDepths,newDistalDepths];

            obj.newPDD    = [newEntrancePDD,newHighResPDD,newDistalPDD];

            %Move to spectra
            nIons = size(obj.ions,2);

            for ionIdx=1:nIons
                eWindow = [obj.phantoms(1).EParam(ionIdx).E(1), obj.phantoms(1).EParam(ionIdx).E(end)];
                nBins = ceil((eWindow(2)-eWindow(1))/obj.downSamplingPhantoms.EParam.Eres);
                %obj.newEnergies(ionIdx,:) = obj.phantoms(1).EParam(1).E;
                obj.newEnergies(ionIdx,:) = linspace(eWindow(1), eWindow(2), nBins);
            end
            
 
             obj.newPhi = cell(1,nIons);
             obj.newedPhi = cell(1,nIons);
             for ionIdx=1:nIons
                 [newInterpDepths,newInterpE] = meshgrid(obj.newDepths,obj.newEnergies(ionIdx,:));
                 [interpDepths,interpEnergies] = meshgrid(obj.depths,obj.phantoms(1).EParam(ionIdx).E);
                 obj.newPhi(ionIdx) = {sparse(interp2(interpDepths,interpEnergies,full(obj.Phi{ionIdx}),newInterpDepths,newInterpE))};
                 obj.newedPhi(ionIdx) = {sparse(interp2(interpDepths,interpEnergies,full(obj.edPhi{ionIdx}),newInterpDepths,newInterpE))};
             end

             % plot profiles
% % %              ionSelect = 6;
% % %              Phi = full(obj.Phi{ionSelect});
% % %              newPhi = full(obj.newPhi{ionSelect});
% % % 
% % %              edPhi = full(obj.edPhi{ionSelect});
% % %              newedPhi = full(obj.newedPhi{ionSelect});
% % %              
% % %              depthSelect = 90; %in mm
% % %              depthSelectIdx = interp1(obj.depths,[1:size(obj.depths,2)],depthSelect, 'nearest');
% % %              newDepthSelectIdx = interp1(obj.newDepths,[1:size(obj.newDepths,2)], depthSelect, 'nearest');
% % %              prof= Phi(:,depthSelectIdx);
% % %              newProf = newPhi(:,newDepthSelectIdx);
% % % 
% % %              profeD= edPhi(:,depthSelectIdx);
% % %              newProfeD = newedPhi(:,newDepthSelectIdx);
% % % 
% % %              eBin = obj.phantoms(1).EParam(1).E;
% % %              newEBin = newEnergies(ionSelect,:);
% % %              figure;
% % %              plot(eBin, prof, '.-');
% % %              hold on;
% % %              plot(newEBin,newProf, '.-');
% % %              legend('old', 'new');
% % %              grid on;
% % %              grid minor;
% % % 
% % %              figure;
% % %              plot(eBin, profeD, '.-');
% % %              hold on;
% % %              plot(newEBin,newProfeD, '.-');
% % %              legend('old', 'new');
% % %              grid on;
% % %              grid minor;

            
        end

       function set.ions(obj, stringValue)
         obj.ions = stringValue;
         nIons = size(stringValue,2);
         for k=1:nIons
            switch stringValue{k}
               case 'protons'
                  Z = 1;
               case 'He'
                  Z = 2;
               case 'Li'
                  Z = 3;
               case 'Be'
                  Z = 4;
               case 'B'
                  Z = 5;
               case 'C'
                  Z = 6;
            end
            obj.ionsZ = [obj.ionsZ, Z];
         end
       end

       function set.phantoms(obj, phantoms)
           obj.phantoms = phantoms;
        for k=1:obj.nPhantoms
             nIons = size(obj.ions,2);
             
             for m=1:nIons
                obj.phantoms(k).EParam(m).EBinWidth = (obj.phantoms(k).EParam(m).EMax-obj.phantoms(k).EParam(m).EMin)/obj.phantoms(k).EParam(m).nEBins;
                obj.phantoms(k).EParam(m).EdBinWidth = (obj.phantoms(k).EParam(m).EdMax-obj.phantoms(k).EParam(m).EdMin)/obj.phantoms(k).EParam(m).nEdBins;
                obj.phantoms(k).EParam(m).E = linspace(obj.phantoms(k).EParam(m).EMin,...
                                                       obj.phantoms(k).EParam(m).EMax - obj.phantoms(k).EParam(m).EBinWidth,...
                                                       obj.phantoms(k).EParam(m).nEBins)...
                                                       +obj.phantoms(k).EParam(m).EBinWidth/2;
     
                obj.phantoms(k).EParam(m).Ed = linspace(obj.phantoms(k).EParam(m).EdMin,...
                                                        obj.phantoms(k).EParam(m).EdMax-obj.phantoms(k).EParam(m).EdBinWidth,...
                                                        obj.phantoms(k).EParam(m).nEdBins)...
                                                        +obj.phantoms(k).EParam(m).EdBinWidth/2;
    
             end
        end

       end
    end
end