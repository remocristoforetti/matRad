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
    end

    methods
        function obj = highResPhantomAnalisis()

        end

        function importRawData(obj)
            depths = [];
            for k=1:obj.nPhantoms
                currDepths = [0:obj.phantoms(k).resolution(1):obj.phantoms(k).dimension(1)-obj.phantoms(k).resolution(1)] + obj.phantoms(k).resolution(1)/2 + obj.phantoms(k).positionDepth(k);
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