classdef DS_MultiIon < handle
   
   properties
      d; %depth in water
      
      Phi;
      rn;
      Beta_Tissue;
      RN;
      LUT = 'D:\matRad_gitHubRemo_MKM_TOPAS\MKM\LUTs_rn25_RN_39_Beta005_500MeV';
      zD; %This is the profile;
      wDir;
      
      EParam;
      Ions;
      
      zD_SingleIons;
 
      IonNumber = [];
      SparsePhi;
      CumulativeEdPhi;
      CumulativePhi;
      PDD;
   end
   methods
      
      function obj = DS_MultiIon()
         
      end
      
      function obj = analyzeRawData(obj, Run)
         matRad_cfg = MatRad_Config.instance();
         nIons = size(obj.EParam,2);

         Phi = [];
         SparsePhi = [];
%          matRad_cfg.dispWarning('Warning! Cutoff on depth inserted');

         CumulativeEdPhi = [];
         CumulativePhi   = [];
         for k=1:nIons
            dataName = strcat(obj.wDir, filesep, 'Ion_', num2str(obj.IonNumber(k)),filesep, 'ZMixOutput_', obj.Ions{k} ,'_Run_', compose('%04d', Run), '.dcm');

            RawD = double(squeeze(dicomread(dataName{1})));
            RawD_infos = dicominfo(dataName{1});
            RawD = RawD*RawD_infos.DoseGridScaling; %->Should scale ion weights relative to each other. Total normazlization to nParicles should not matter
            RawD = permute(RawD, [2,3,1]);
            %Phi = [Phi, {RawD(2:end,2:end,:)}];

           R = RawD(2:end,2:end,:);
           Ed = obj.EParam(k).Ed;
           Edp = repmat(Ed,1,1,size(R,3));
           CumulativePhi = [CumulativePhi, {sparse(squeeze(sum(R, [1])))}];
           CumulativeEdPhi = [CumulativeEdPhi, {sparse(squeeze(pagemtimes(Edp,R)))}];
           %SparsePhi = [SparsePhi, {sparse(R(:))}];
              
         end
         %obj.SparsePhi = SparsePhi;
         obj.CumulativeEdPhi = CumulativeEdPhi;
         obj.CumulativePhi   = CumulativePhi;
         %obj.Phi = Phi;

      end
      
      function PDD = readDoseProfile(obj,energy,nPixels)

         fileName = dir(strcat(obj.wDir, '*_physicalDose_Run_', string(compose('%04d', energy)), '.dcm'));
         rawProfile = double(squeeze(dicomread(strcat(fileName.folder,filesep,fileName.name))));
         infos = dicominfo(strcat(fileName.folder,filesep,fileName.name));
         PDD = rawProfile.*infos.DoseGridScaling;
         
         PDD = permute(PDD, [2,3,1]);
         integIntervall_1  = [ceil(size(PDD,1)/2)-nPixels:ceil(size(PDD,1)/2)+nPixels];
         integIntervall_2  = [ceil(size(PDD,2)/2)-nPixels:ceil(size(PDD,2)/2)+nPixels];
         
         obj.PDD = squeeze(sum(PDD(integIntervall_1,integIntervall_2,:), [1,2]));
      end
      
       function obj = computeZD(obj)         
         nIons = size(obj.EParam,2);

         zD = zeros(size(obj.d,2),1);
         totEdep = zeros(size(obj.d,2),1);
         for k = 1:nIons
            lut = load(strcat(obj.LUT, filesep, 'LUT_Z_', num2str(k), '.mat'));
            ELUT = lut.E;
            zDLUT = lut.zD;
            zDk_ion = interp1(ELUT,zDLUT,obj.EParam(k).E);
            zDk_ion(isnan(zDk_ion)) = 0;
            %Maybe for loop not needed
            %zD = [];
            for m=1:size(obj.d,2)
               totEdep(m) = totEdep(m) + sum(obj.EParam(k).Ed*obj.Phi{k}(:,:,m));
               zD(m) = zD(m) + zDk_ion*(obj.EParam(k).Ed*obj.Phi{k}(:,:,m))';
            end
         end
         obj.zD = zD./totEdep;
       end
 
       function obj = computeZD_SingleIons(obj)
         nIons = size(obj.EParam,2);
         obj.zD_SingleIons = [];
         
         for k = 1:nIons

            lut = load(strcat(obj.LUT, filesep, 'LUT_Z_', num2str(k), '.mat'));
            ELUT = lut.E;
            zDLUT = lut.zD;
            zDk_ion = interp1(ELUT,zDLUT,obj.EParam(k).E);
         
            %Maybe for loop not needed
            %zD = [];
            for m=1:size(obj.d,2)
               totEdep = sum(obj.EParam(k).Ed*obj.Phi{k}(:,:,m));
               if totEdep > 0
                  zD(m) = zDk_ion*(obj.EParam(k).Ed*obj.Phi{k}(:,:,m))'/totEdep;
               else
                  zD(m) = 0; 
               end
            end
            obj.zD_SingleIons = [obj.zD_SingleIons {zD}];
         end
      end
      
      function computeBinning(obj)
         nIons = size(obj.EParam,2);
         
         for k=1:nIons
            obj.EParam(k).EBinWidth = (obj.EParam(k).EMax-obj.EParam(k).EMin)/obj.EParam(k).nEBins;
            obj.EParam(k).EdBinWidth = (obj.EParam(k).EdMax-obj.EParam(k).EdMin)/obj.EParam(k).nEdBins;
            obj.EParam(k).E = linspace(obj.EParam(k).EMin,obj.EParam(k).EMax,obj.EParam(k).nEBins)+obj.EParam(k).EBinWidth/2;
 
            obj.EParam(k).Ed = linspace(obj.EParam(k).EdMin,obj.EParam(k).EdMax,obj.EParam(k).nEdBins)+obj.EParam(k).EdBinWidth/2;

         end
         
      end

      function set.Ions(obj, stringValue)
         obj.Ions = stringValue;
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
            obj.IonNumber = [obj.IonNumber, Z];
         end
      end
   end
end