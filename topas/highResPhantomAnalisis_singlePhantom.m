classdef highResPhantomAnalisis_singlePhantom < handle
    properties
        nphantom;
        wDir;
        readPDD;
        
        readDSSpectra;
        readSpectra;
        readSPeDSpectra;
        readProtonLET;
        readSPdose;
        readEventsCounter
        depths;
        phantom;
        PDD;

        Phi;

        edPhi;
        fPhi;
        ProtonLET;
        DSPhi;
        run;
        integrationSurface;
        ions;
        ionsZ;
        downSamplingPhantom;
        newDepths;

        newPDD;
        newPhi;
        newedPhi;
        newfPhi
        newEnergies;
        newDosePhi;
        EParam;
        spectraLET;
        spectraLETdS;
        spectraLETt;
        dosePhi;
        averageQuantities;
    end

    methods
        function obj = highResPhantomAnalisis_singlePhantom()

        end

        function importRawData(obj)
  
            if obj.readPDD
                scorersDir = [obj.wDir,filesep, 'Output',filesep, 'PDD', filesep, obj.phantom{1,5}.name];
                dataFilename = strcat(scorersDir, filesep, 'physicalDose_Run_', compose('%04d', obj.run), '.dcm');
                data_dicom = double(squeeze(dicomread(dataFilename{1})));
                data_dicom_info = dicominfo(dataFilename{1});
                data_dicom_permute = permute(data_dicom, [2,3,1]);
                data_dicom = squeeze(sum(data_dicom_permute(obj.integrationSurface, obj.integrationSurface, :), [1,2])).*data_dicom_info.DoseGridScaling;
                obj.PDD = data_dicom;
            end

% % %             if obj.readDSSpectra %still to be reviewed, not used
% % %                 scorersDir = [obj.wDir,filesep, 'Output',filesep, 'DS', filesep, 'Phantom'];
% % %                 data = [];
% % % 
% % %                 for m = 1:size(obj.ionsZ,2)
% % %                     Phi = [];
% % %                     edPhi = [];
% % %                     for k = 1:obj.nphantom
% % %                         dataFilename = strcat(scorersDir, num2str(k), filesep, 'Ion_', num2str(obj.ionsZ(m)),filesep, 'S_',obj.ions{m}, '_Run_', compose('%04d', obj.run), '.dcm');
% % %                         data_dicom = double(squeeze(dicomread(dataFilename{1})));
% % %                         data_dicom_info = dicominfo(dataFilename{1});
% % %                         data_dicom_permute = permute(data_dicom, [2,3,1]);
% % %                         data_dicom = data_dicom_permute(2:end,2:end,:).*data_dicom_info.DoseGridScaling;%./obj.phantom(k).resolution(1); %This is corrected for phantom depth resolution. In the end, computation of alpha/beta does not depend on normalization
% % %                     
% % %                         Ed = obj.phantom(k).EParam(m).Ed;
% % %                         Ed = repmat(Ed,1,1,size(data_dicom,3));
% % %                         
% % %                         Phi = [Phi, sparse(squeeze(sum(data_dicom, [1])))];
% % %                         edPhi = [edPhi, sparse(squeeze(pagemtimes(Ed,data_dicom)))];
% % %                     
% % %                     end
% % %                     obj.DSPhi{m} = Phi;
% % % 
% % %                     obj.edPhi{m} = edPhi;
% % %                 end
% % % 
% % %             end

            %These are the energy spectra, count for number of ions, not
            %number of events
            if obj.readSpectra
                scorersDir = [obj.wDir,filesep, 'Output',filesep, 'SP', filesep, obj.phantom{1,5}.name];

                for m = 1:size(obj.ionsZ,2)
                    dataFilename = strcat(scorersDir,filesep, 'Ion_', num2str(obj.ionsZ(m)),filesep, 'SP_',obj.ions{m}, '_Run_', compose('%04d', obj.run), '.csv');
                    data_csv = csvread(dataFilename{1},11);
                    %This is fluence, (1/mm^2), su multiply by phantom
                    %lateral surface (scoring is performed over whole phantom surface)
                    data_csv = data_csv(:,2:end-2)'.*prod([obj.phantom{1,5}.dimension.x,obj.phantom{1,5}.dimension.z]);
                    obj.fPhi{m} = sparse(squeeze(data_csv));
                end

            end

            if obj.readEventsCounter

                scorersDir = [obj.wDir,filesep, 'Output',filesep, 'Phi', filesep, obj.phantom{1,5}.name];
                for m = 1:size(obj.ionsZ,2)
                    dataFilename = strcat(scorersDir, filesep, 'Ion_', num2str(obj.ionsZ(m)),filesep, 'Phi_',obj.ions{m}, '_Run_', compose('%04d', obj.run), '.dcm');
                    data_dicom = double(squeeze(dicomread(dataFilename{1})));
                    data_dicom_info = dicominfo(dataFilename{1});
                    data_dicom = data_dicom*data_dicom_info.DoseGridScaling;%./obj.phantom{1,5}.resolution.y;
                    obj.Phi{m} = sparse(squeeze(data_dicom(:,2:end)'));
                end

            end



            if obj.readSPeDSpectra

                scorersDir = [obj.wDir,filesep, 'Output',filesep, 'SPeD', filesep, obj.phantom{1,5}.name];
                for m = 1:size(obj.ionsZ,2)
                    dataFilename = strcat(scorersDir, filesep, 'Ion_', num2str(obj.ionsZ(m)),filesep, 'SPeD_',obj.ions{m}, '_Run_', compose('%04d', obj.run), '.dcm');
                    data_dicom = double(squeeze(dicomread(dataFilename{1})));
                    data_dicom_info = dicominfo(dataFilename{1});
                    data_dicom = data_dicom*data_dicom_info.DoseGridScaling;%./obj.phantom{1,5}.resolution.y;
                    

                    if obj.averageQuantities
                        obj.edPhi{m} = sparse(squeeze(data_dicom(:,2:end)'));
                    else
                        obj.edPhi{m} = sparse(squeeze(data_dicom(:,2:end)').*obj.Phi{m});
                    end
                end

            end

             if obj.readProtonLET %To be reviewed, not used
                scorersDir = [obj.wDir,filesep, 'Output',filesep, 'ProtonLET', filesep, 'Phantom'];
                data = [];
                LET = [];

                for k = 1:obj.nphantom

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

             %%%compute LET depth profle
            

           if obj.readSPdose
                scorersDir = [obj.wDir,filesep, 'Output',filesep, 'SPdose', filesep, obj.phantom{1,5}.name];
                for m = 1:size(obj.ionsZ,2)
                    dataFilename = strcat(scorersDir, filesep, 'Ion_', num2str(obj.ionsZ(m)),filesep, 'SPdose_',obj.ions{m}, '_Run_', compose('%04d', obj.run), '.dcm');
                    data_dicom = double(squeeze(dicomread(dataFilename{1})));
                    data_dicom_info = dicominfo(dataFilename{1});

                    data_dicom = data_dicom*data_dicom_info.DoseGridScaling;
                    if obj.averageQuantities
                        obj.dosePhi{m} = sparse(squeeze(data_dicom(:,2:end)'));
                    else
                        obj.dosePhi{m} = sparse(squeeze(data_dicom(:,2:end)').*obj.Phi{m});
                    end
                end

           end
           if ~isempty(obj.dosePhi)

               obj.spectraLET = cell(1,1);
               obj.spectraLETt = cell(1,1);
               fullLETDepth = 0;
               fullLETDeptht = 0;
               normDenominator = 0;
               normDenominatort = 0;
                for ionIdx =1:size(obj.ionsZ,2)
                    letEnergy = arrayfun(@(e) obj.computeLET(e,ionIdx), obj.phantom{1,4}(ionIdx).E);
                    
                    
                    % LETd
                    fullLETDepth = fullLETDepth+ letEnergy*(full(obj.dosePhi{ionIdx}));%.*full(obj.Phi{ionIdx}));
                    fullLETDeptht = fullLETDeptht + letEnergy*full(obj.fPhi{ionIdx});
                
                    normDenominator = normDenominator + sum(full(obj.dosePhi{ionIdx}));%.*full(obj.Phi{ionIdx}),[1]);
                    normDenominatort = normDenominatort + sum(full(obj.fPhi{ionIdx}),[1]);
                
                end

                depthsToCompute = normDenominator >0;
                depthsToComputet = normDenominatort >0;
                
                obj.spectraLET{1}(depthsToCompute) = fullLETDepth(depthsToCompute)./normDenominator(depthsToCompute);
                obj.spectraLETt{1}(depthsToComputet) = fullLETDeptht(depthsToComputet)./normDenominatort(depthsToComputet);
                
                

%                obj.spectraLET = cell(1,size(obj.ionsZ,2));
%                obj.spectraLETt = cell(1,size(obj.ionsZ,2));
%                  for ionIdx =1:size(obj.ionsZ,2)
%                     letEnergy = arrayfun(@(e) obj.computeLET(e,ionIdx), obj.phantom{1,4}(ionIdx).E);
%                     obj.spectraLET{ionIdx} = zeros(1,size(obj.depths,2));
%                     obj.spectraLETt{ionIdx} = zeros(1,size(obj.depths,2));
%                     
%                     % LETd
%                     fullLETDepth = letEnergy*(full(obj.dosePhi{ionIdx}).*full(obj.Phi{ionIdx}));
%                     normDenominator = sum(full(obj.dosePhi{ionIdx}).*full(obj.Phi{ionIdx}),[1]);
%                     depthsToCompute = normDenominator > 0;
% 
%                     obj.spectraLET{ionIdx}(depthsToCompute) = fullLETDepth(depthsToCompute)./normDenominator(depthsToCompute);
%  %                   obj.spectraLET{ionIdx} = fullLETDepth;
% 
%                     % LETt
% 
%                     fullLETDeptht = letEnergy*full(obj.Phi{ionIdx});
%                     normDenominatort = sum(full(obj.Phi{ionIdx}),[1]);
% 
%                     depthsToComputet = normDenominatort > 0;
% 
%                     obj.spectraLETt{ionIdx}(depthsToComputet) = fullLETDeptht(depthsToComputet)./normDenominatort(depthsToComputet);
%                     %obj.spectraLETt{ionIdx} = fullLETDepth;
%                  end
            end
            obj.rawDataDownsampling();
        end





        function rawDataDownsampling(obj)
   
            res = obj.downSamplingPhantom.resolutions;
      

            highResWindowWidth = obj.downSamplingPhantom.highResWindowWidth;

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
            entranceWindowWidth = [obj.depths(1),highResWindow(1)];
            nBins = ceil(entranceWindowWidth(2)/res(1));

            newEntranceDepths = linspace(entranceWindowWidth(1),entranceWindowWidth(2),nBins+1);
            newEntranceDepths = newEntranceDepths(1:end-1);
            newEntrancePDD = interp1([obj.depths], [obj.PDD], newEntranceDepths,'spline');

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
            obj.newPhi = cell(1,nIons);
            obj.newedPhi = cell(1,nIons);
            obj.newfPhi = cell(1,nIons);
            obj.newDosePhi = cell(1,nIons);
            for ionIdx=1:nIons
               nBins = obj.downSamplingPhantom.EParam(ionIdx).nEnergies;
               eLine = [];
               if obj.downSamplingPhantom.EParam(ionIdx).eLine(2) >= 0
 
                   eLine = obj.downSamplingPhantom.EParam(ionIdx).eLine.*ones(size(obj.newDepths,2),1); %then this can change and be provided from external as a function of depth

               else
                %This makes sense only if distribution is more or less a
                %gaussian

                %Get rough index of the maximum of the distributuion
                [maxValue,maxIdx] = max(full(obj.Phi{ionIdx}));

                %Get energy of the max
                maxE = obj.phantom{1,4}(ionIdx).E(maxIdx);

                %Setup an arbitrary window in MeV/u
                eWindowWidth = 200; %MeV/u;
                eLineOffset = [];
                eLineOffset= maxE - eWindowWidth/2;
                eLineOffset(eLineOffset<0) = obj.phantom{1,4}(ionIdx).E(1);

                eLine(:,1) = interp1(obj.depths,eLineOffset,obj.newDepths');
                maxEbinnig = interp1(obj.depths,maxE + eWindowWidth/2,obj.newDepths);
                eLine(:,2) = (maxEbinnig - eLine(:,1)')./nBins;

               end
               
               obj.EParam(ionIdx).nEnergies = nBins;
               obj.EParam(ionIdx).eLine     = eLine;
               newInterpE = (eLine(:,2)*[0:nBins-1] + eLine(:,1))';
               
                newInterpDepths = repmat(obj.newDepths,nBins,1);

               [interpDepths,interpEnergies] = meshgrid(obj.depths,obj.phantom{1,4}(ionIdx).E);

               obj.newPhi(ionIdx) = {sparse(interp2(interpDepths,interpEnergies,full(obj.Phi{ionIdx}),newInterpDepths,newInterpE))};
               obj.newedPhi(ionIdx) = {sparse(interp2(interpDepths,interpEnergies,full(obj.edPhi{ionIdx}),newInterpDepths,newInterpE))};
               obj.newfPhi(ionIdx) = {sparse(interp2(interpDepths,interpEnergies,full(obj.fPhi{ionIdx}),newInterpDepths,newInterpE))};
               obj.newDosePhi(ionIdx) = {sparse(interp2(interpDepths,interpEnergies,full(obj.dosePhi{ionIdx}),newInterpDepths,newInterpE))};
            end
            
            %%downsample LET profiles

%             for ionIdx = 1:nIons
%                 obj.spectraLETdS{ionIdx} = interp1(obj.depths, obj.spectraLET{ionIdx}, obj.newDepths, 'spline');
%             end

              obj.spectraLETdS = interp1(obj.depths, obj.spectraLET{1}, obj.newDepths, 'spline');

%             ionSelect = 6;
%             eIdx = 1;
%             Phi = full(obj.Phi{ionSelect});
%             newPhi = full(obj.newPhi{ionSelect});
%             
%             edPhi = full(obj.edPhi{ionSelect});
%             newedPhi = full(obj.newedPhi{ionSelect});
%             
%             depthSelect = 10; %in mm
%             
%             depthSelectIdx = interp1(obj.depths,[1:size(obj.depths,2)],depthSelect, 'nearest');
%             newDepthSelectIdx = interp1(obj.newDepths,[1:size(obj.newDepths,2)], depthSelect, 'nearest');
%             prof= Phi(:,depthSelectIdx);
%             newProf = newPhi(:,newDepthSelectIdx);
%             
%             profeD= edPhi(:,depthSelectIdx);
%             newProfeD = newedPhi(:,newDepthSelectIdx);
%             EBinWidth = (obj.phantom{1,4}(ionSelect).EMax-obj.phantom{1,4}(ionSelect).EMin)/obj.phantom{1,4}(ionSelect).nEBins;
%             eBin = obj.phantom{1,4}(ionIdx).E
%             
%             newEBin = newInterpE(:,newDepthSelectIdx);
%             figure;
%             
%             plot(eBin, prof, '.-');
%             hold on;
%             plot(newEBin,newProf, '.-');
%             legend('old', 'new');
%             grid on;
%             grid minor;


            
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

       function set.phantom(obj, phantom)
        obj.phantom = phantom;
        nIons = size(obj.ions,2);
        for k=1:nIons
            obj.phantom{1,4}(k).EBinWidth = (obj.phantom{1,4}(k).EMax-obj.phantom{1,4}(k).EMin)/obj.phantom{1,4}(k).nEBins;
            obj.phantom{1,4}(k).E = linspace(obj.phantom{1,4}(k).EMin,...
                                                   obj.phantom{1,4}(k).EMax - obj.phantom{1,4}(k).EBinWidth,...
                                                   obj.phantom{1,4}(k).nEBins)...
                                                   +obj.phantom{1,4}(k).EBinWidth/2;
         end

         nBins = ceil(phantom{1,5}.dimension.y/phantom{1,5}.resolution.y);
         obj.depths = linspace(0,phantom{1,5}.dimension.y-phantom{1,5}.resolution.y, nBins) + phantom{1,5}.resolution.y/2;
       end

         function let = computeLET(obj,E,Zion)
         
         %Constants
         I = 75*10^(-6);
         A_1 = 4.382*10^(-25); %(8.4018*10^(-54))/(4*pi*0.51); n*e^4/4*pi*mec^2*e0^2; %https://www.nature.com/articles/s41598-017-10554-0
         A_2 = 0.511; %mec^2
         A_3 = 3.893*10^(22);
         AMU2MEV = 931.494;
         %Rest energy 
         switch Zion
                case 1 %p
                    A = 1;
                case 2 %He
                    A = 4;
                case 3  %Li
                    A = 7;
                case 4  %Be
                    A = 9;
             case 5  %B %%%% Survival uses A = 10 for B %%%%
                    A = 11;
                case 6  %C
                    A = 12;
         end
    
         Erest = A*AMU2MEV;
         E = A*E;
         Beta_Ion = sqrt(1 - 1./((E./Erest +1).^2));
         
         LOG = log((2*A_2/I)*((Beta_Ion.^2)./(1-Beta_Ion.^2)));
         let = A_1*((Zion./Beta_Ion).^2).*(LOG - Beta_Ion.^2)*A_3;
         end

         function eMax = maxEnergy(obj, ion)
            eLine = obj.downSamplingPhantom.EParam(ion).eLine;
            nEnergies = obj.downSamplingPhantom.EParam(ion).nEnergies;
            eMax = eLine(2)*(nEnergies-1) + eLine(1);
         end
    end
end