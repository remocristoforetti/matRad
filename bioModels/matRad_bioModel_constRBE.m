classdef matRad_bioModel_constRBE < matRad_BiologicalModel
      properties (SetAccess = private, Hidden)%constant
        constRBE_protons = 1.1;
        constRBE_photons = 1;
        default_AlphaX = 0.1;
        default_BetaX = 0.05;
        
        default_RBE_protons = 1.1;
        default_RBE_photons = 1;
        
        AvailableRadiationModalities = {'protons'};
        RequiredBaseData = {};
        RBE;
      end

      %% Methods
   methods
       %Constructor
       function obj = matRad_bioModel_constRBE(radiationMode)
         if ~ischar(radiationMode)
           matRad_dispToConsole(['Something wrong with bioModel inputs'],[],'warning');
         end
         
         obj@matRad_BiologicalModel(radiationMode);
         obj.model = 'constRBE';

         switch obj.radiationMode
             % case {'photons'}
             %       obj.RBE = obj.default_RBE_photons;
             case {'protons'}
                   obj.RBE = obj.default_RBE_protons;
         end
       end
 
       function str = calcTissueParameters(obj,cst,numVoxels,stf,ctScen)
           matRad_cfg = MatRad_Config.instance();
           str = struct('alphaX', [], ...
                        'betaX', [], ...
                        'tissueABratio', [], ...
                        'RBE', [], ...
                        'tissueClass', []);
           
           TissueParam = {cst{:,5}};
           tissueClass = zeros(numVoxels,1);
           for k = 1:size(TissueParam,2)
               
              if ~isfield(TissueParam{k}, 'alphaX')
                  str.alphaX = [str.alphaX, obj.default_AlphaX];
                  matRad_cfg.dispWarning(['Warning! Tissue parameter AlphaX for tissue %i not available, set do default'], cst{k,1});
              else
                  str.alphaX = [str.alphaX, TissueParam{k}.alphaX];
              end
               
              if ~isfield(TissueParam{k}, 'betaX')
                  str.betaX = [str.betaX, obj.default_BetaX];
                  matRad_cfg.dispWarning(['Warning! Tissue parameter BetaX for tissue %i not available, set do default '],  cst{k,1});
              else
                  str.betaX = [str.betaX, TissueParam{k}.betaX];
              end
              
              if ~isfield(TissueParam{k}, 'RBE')
                  switch obj.radiationMode

                     case {'protons'}
                        str.RBE = [str.RBE, obj.default_RBE_protons];
                  end
                  matRad_cfg.dispWarning(['Warning! Tissue parameter RBE for tissue not available, set to default'], cst{k,1});
              else
                  str.RBE = [str.RBE, TissueParam{k}.RBE];
              end

              
              str.tissueABratio(k) = str.alphaX(k)./str.betaX(k);
              tissueClass(cst{k,4}{1},1) = k;
           end
           str.tissueClass = tissueClass(tissueClass>0);
       end
       
       function [bixelAlpha,bixelBeta] = calcLQParameter(obj,vRadDepths,baseDataEntry,TissueParam,ix)
            
            bixelAlpha = NaN*ones(numel(vRadDepths),1);
            bixelBeta  = NaN*ones(numel(vRadDepths),1);
            
            %This is done so that if proj_variableRBE is called with this
            %model, it gives back RBExD anyway
            RBE = TissueParam.RBE(TissueParam.tissueClass(ix))';
            alphaX = TissueParam.alphaX(TissueParam.tissueClass(ix))';
            betaX = TissueParam.betaX(TissueParam.tissueClass(ix))';
            
            bixelAlpha = RBE.*alphaX;
            bixelBeta  = (RBE.^2).*betaX;
                    
       end
   end
   
end