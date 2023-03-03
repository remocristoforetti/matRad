classdef matRad_BioModel_WED < matRad_Bio_Model
%% Properites
    properties (SetAccess = private)%constant
        p0_WED   = 1; % https://www.ncbi.nlm.nih.gov/pubmed/22909391
        p1_WED   = 0.434;
        p2_WED   = 1;
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
       %Constructor
       function obj = matRad_BioModel_WED(radiationMode)
         if ~ischar(radiationMode)
           matRad_cfg.dispWarning(['Something wrong with bioModel inputs']);
          end
          
          obj@matRad_Bio_Model();
          obj.model = 'WED';
          obj.radiationMode = radiationMode;

       end
       
       function str = calcTissueParameter(obj,cst,numVoxels,ctScen)
           
           str = struct('alphaX', [], ...
                        'betaX', [], ...
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
       
       function [bixelAlpha,bixelBeta] = calcLQParameter(obj,vRadDepths,baseDataEntry,TissueParam, ix)
            
            bixelAlpha = NaN*ones(numel(vRadDepths),1);
            bixelBeta  = NaN*ones(numel(vRadDepths),1);
            
            % range shift
            depths = baseDataEntry.depths + baseDataEntry.offset;
 
            bixelLET = matRad_interp1(depths,baseDataEntry.LET,vRadDepths);
            bixelLET(isnan(bixelLET)) = 0;

            vAlpha_x = TissueParam.alphaX(TissueParam.tissueClass(ix))';
            vBeta_x = TissueParam.betaX(TissueParam.tissueClass(ix))';
            vABratio = TissueParam.tissueABratio(TissueParam.tissueClass(ix))';

            RBEmax     = obj.p0_WED + ((obj.p1_WED * bixelLET )./ vABratio);
            RBEmin     = obj.p2_WED;
            bixelAlpha = RBEmax    .* vAlpha_x;
            bixelBeta  = RBEmin.^2 .* vBeta_x;
                    
       end
   end
end