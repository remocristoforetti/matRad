classdef matRad_BioModel_CAR < matRad_Bio_Model
%% Properties
   properties (SetAccess = private)
        p0_CAR   = 0.843; % http://www.tandfonline.com/doi/abs/10.1080/09553000601087176?journalCode=irab20
        p1_CAR   = 0.154; % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4869317/pdf/13014_2016_Article_642.pdf
        p2_CAR   = 2.686;
        p3_CAR   = 1.09;
        p4_CAR   = 0.006;
        p5_CAR   = 2.686;
      
      
      
        default_AlphaX = 0.1;
        default_BetaX = 0.5;
        
        AvailableradiationModalities = {'protons'};
        AvailableQuantitiesForOpt = {'physicalDose','RBExD','effect'};
        RequiredBaseData = {'depths','offset','LET'};
        NumModelParameters = 2;
        ModelParameters = {'alphaX','betaX'};
      
   end
%% Methods
   methods
      function obj = matRad_BioModel_CAR(radiationMode)

         matRad_cfg = MatRad_Config.instance();
         if ~ischar(radiationMode)
              matRad_cfg.dispWarning(['Something wrong with bioModel inputs \n']);
         end
         obj@matRad_Bio_Model();
         obj.model = 'CAR';
         obj.radiationMode = radiationMode;
      end


      function str = calcTissueParameter(obj,cst,numVoxels,ctScen)
         matRad_cfg = MatRad_Config.instance();  
         str = struct('alphaX', [], ...
                      'betaX', [], ...
                      'tissueABratio', [], ...
                      'tissueClass', []);
         TissueParam = {cst{:,5}};
          
         tissueClass = zeros(numVoxels,1);
         for k = 1:size(TissueParam,2)
              
            if ~isfield(TissueParam{k}, 'alphaX')
                str.alphaX = [str.alphaX, obj.default_AlphaX];
                matRad_cfg.dispWarning(['Warning! Tissue parameter AlphaX for tissue ' cst{k,1} 'not available, set do default \n']);
            else
                str.alphaX = [str.alphaX, TissueParam{k}.alphaX];
            end
              
            if ~isfield(TissueParam{k}, 'betaX')
                str.betaX = [str.betaX, obj.default_BetaX];
                matRad_cfg.dispWarning(['Warning! Tissue parameter BetaX for tissue ' cst{k,1} 'not available, set do default \n']);
            else
                str.betaX = [str.betaX, TissueParam{k}.betaX];
            end
              
            str.tissueABratio(k) = str.alphaX(k)./str.betaX(k);
            tissueClass(cst{k,4}{1},1) = k; %What if there is an overlap between structures?
         end
         str.tissueClass = tissueClass(tissueClass>0);

      end
      
      function [bixelAlpha,bixelBeta] = calcLQParameter(obj,vRadDepths,baseDataEntry,TissueParam,ix)
            
         bixelAlpha = NaN*ones(numel(vRadDepths),1);
         bixelBeta  = NaN*ones(numel(vRadDepths),1);
            
            % range shift
         depths = baseDataEntry.depths + baseDataEntry.offset;
 
         bixelLET = matRad_interp1(depths,baseDataEntry.LET,vRadDepths);
         bixelLET(isnan(bixelLET)) = 0;
             
         vAlpha_x = TissueParam.alphaX(TissueParam.tissueClass(ix))';
         vBeta_x = TissueParam.betaX(TissueParam.tissueClass(ix))';
         vABratio = TissueParam.tissueABratio(TissueParam.tissueClass(ix))';
         
         gamma      = obj.p2_CAR./vABratio;
         RBEmax     = obj.p0_CAR + obj.p1_CAR.*gamma.*bixelLET;
         RBEmin     = obj.p3_CAR + obj.p4_CAR.*gamma.*bixelLET;
             
             
         bixelAlpha = RBEmax    .* vAlpha_x;
         bixelBeta  = RBEmin.^2 .* vBeta_x;
      end
   end
end