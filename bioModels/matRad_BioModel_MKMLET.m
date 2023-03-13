classdef matRad_BioModel_MKMLET < matRad_Bio_Model

   properties (SetAccess = private, Hidden)%constant
       
       %As taken from TOPAS implementation of MKMLET model (https://github.com/topasmc/extensions)
       default_rD = 0.35;%0.52 %+ 0.5*0.52;        %0.52 um, domain diameter only valid for V79 Cells! (TOPAS)
       default_Alpha0 = 0.1%Gy^(-1), only valid for V79 Cells! (TOPAS)
       default_AlphaX = 0.1;
       default_BetaX = 0.05;
       default_Beta0 = 0.05;
       Gamma0 = 0.229   %Constant value Following Hawkins 1997 (https://doi.org/10.1118/1.598307)
       
       AvailableradiationModalities = {'protons'};
       AvailableQuantitiesForOpt = {'physicalDose','RBExD','effect'};
       RequiredBaseData = {'depths','offset','LET'};
       NumModelParameters = 5;
       ModelParameters = {'alpha0', 'beta0','rD','alphaX','betaX'};
   end

      %% Methods
   methods

       %Constructor
       function obj = matRad_BioModel_MKMLET(radiationMode)
          if ~ischar(radiationMode)
           matRad_cfg.dispWarning(['Something wrong with bioModel inputs \n']);
          end

          obj@matRad_Bio_Model();
          obj.model = 'MKMLET';
          obj.radiationMode = radiationMode;
       end
       
       function str = calcTissueParameters(obj,cst,numVoxels,stf,ctScen)
           matRad_cfg = MatRad_Config.instance();
           str = struct('alpha0', [], ...
                        'beta0', [], ...
                        'rD', [], ...
                        'alphaX', [], ...
                        'betaX', [], ...
                        'tissueClass', []);
           
           obj.NumModelParameters = sum(~strcmp(fieldnames(str), 'tissueClass'));
           
           TissueParam = {cst{:,5}};
           tissueClass = zeros(numVoxels,1);
           for k = 1:size(TissueParam,2)
               
              if ~isfield(TissueParam{k}, 'alphaX')
                  str.alphaX = [str.alphaX, obj.default_AlphaX];
                  matRad_cfg.dispWarning(['Warning! Tissue parameter AlphaX for tissue %i  not available, set do default'],cst{k,1});
              else
                  str.alphaX = [str.alphaX, TissueParam{k}.alphaX];
              end
               
              if ~isfield(TissueParam{k}, 'betaX')
                  str.betaX = [str.betaX, obj.default_BetaX];
                  matRad_cfg.dispWarning(['Warning! Tissue parameter BetaX for tissue %i not available, set do default'],cst{k,1});
              else
                  str.betaX = [str.betaX, TissueParam{k}.betaX];
              end

              if ~isfield(TissueParam{k}, 'alpha0')
                  str.alpha0 = [str.alpha0, obj.default_Alpha0];
                  matRad_cfg.dispWarning(['Warning! Tissue parameter Alpha0 for tissue %i not available, set do default'], cst{k,1});
              else
                  str.alpha0 = [str.alpha0, TissueParam{k}.alpha0];
              end
 
              if ~isfield(TissueParam{k}, 'beta0')
                  str.beta0 = [str.beta0, obj.default_Beta0];
                  matRad_cfg.dispWarning(['Warning! Tissue parameter Beta0 for tissue %i not available, set do default'], cst{k,1});
              else
                  str.beta0 = [str.beta0, TissueParam{k}.beta0];
              end

              if ~isfield(TissueParam{k}, 'rD')
                  str.rD = [str.rD, obj.default_rD];
                  matRad_cfg.dispWarning(['Warning! Tissue parameter rD for tissue %i not available, set do default'], cst{k,1});
              else
                  str.rD = [str.rD, TissueParam{k}.rD];
              end
              
              tissueClass(cst{k,4}{1},1) = k; %What if there is an overlap between structures?
           end
           str.tissueClass = tissueClass(tissueClass>0);

       end
       
       function [bixelAlpha,bixelBeta] = calcLQParameter(obj,vRadDepths,baseDataEntry,TissueParam,ix)%calcLQParameter(obj,vRadDepths,baseDataEntry,mTissueClass,tissueParameters),vAlpha_x,vBeta_x,vABratio);
            
%             if ~isfield(tissueParameters,'alpha0')
%                 alpha0 = tissueParameters.alphaX
           
            bixelAlpha = NaN*ones(numel(vRadDepths),1);
            bixelBeta  = NaN*ones(numel(vRadDepths),1);
            
            % range shift
            depths = baseDataEntry.depths + baseDataEntry.offset;
 
            bixelLET = matRad_interp1(depths,baseDataEntry.LET,vRadDepths);
            bixelLET(isnan(bixelLET)) = 0;
            
            rD = TissueParam.rD(TissueParam.tissueClass(ix))';
            Alpha0 = TissueParam.alpha0(TissueParam.tissueClass(ix))';
            Beta0 = TissueParam.beta0(TissueParam.tissueClass(ix))';
            
            gamma = (obj.Gamma0./(rD.^2)).*bixelLET;
            
            bixelAlpha = Alpha0 + gamma.*Beta0;
            bixelBeta = Beta0;
       end
   end
end