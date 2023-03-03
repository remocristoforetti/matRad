classdef matRad_BioModel_LEM < matRad_Bio_Model

    properties (SetAccess = private)%constant
  
        AvailableradiationModalities = {'carbon', 'helium'};
        AvailableQuantitiesForOpt = {'physicalDose','RBExD','effect'};
        RequiredBaseData = {'depths','offset','alpha', 'beta', 'alphaX', 'betaX'};
        default_AlphaX = 0.1;
        default_BetaX = 0.05;
   end

      %% Methods

      methods
       %Constructor
       function obj = matRad_BioModel_LEM(radiationMode)
          if ~ischar(radiationMode)
           matRad_cfg.dispWarning(['Something wrong with bioModel inputs']);
          end
          
          obj@matRad_Bio_Model();
          obj.model = 'LEM';
          obj.radiationMode = radiationMode;
       end
       

       function str = calcTissueParameters(obj,cst,numVoxels,stf,ctScen)
           matRad_cfg =  MatRad_Config.instance();

           str = struct('tissueClass', []);
           
           TissueParam = {cst{:,5}};
           tissueClass = zeros(numVoxels,1);
           
           try
              load([matRad_cfg.matRadRoot filesep 'basedata' filesep strcat(obj.radiationMode, '_', obj.baseData) '.mat']);
           catch
              matRad_cfg.dispError(['Could not find the following machine file: ' obj.machine '\n']);
           end

           if exist('machine', 'var')
              for k = 1:size(TissueParam,2)

%                  if ~isfield(TissueParam{k}, 'alphaX') || (isfield(TissueParam{k}, 'alphaX') && ~any(ismember(machine.data(1).alphaX,TissueParam{k}.alphaX)))
%                      %Is matlab evaluating the whole condition anyway? Does
%                      %it stop if first condition is met?
%                      matRad_cfg.dispWarning(['Warning! Tissue parameter AlphaX for tissue ' cst{k,1} 'not available, set do default \n']);
%                      if ismember(obj.default_AlphaX, machine.data(1).alphaX)
%                          str.alphaX = [str.alphaX, obj.default_AlphaX];
%                      else
%                          str.alphaX = [str.alphaX, machine.data(1).alphaX(1)];
%                          matRad_cfg.dispWarning(['Warning! Tissue parameter default_AlphaX for tissue ' cst{k,1} 'not compatible with basedata, setting the available one \n']);
%                      end
%                  else
%                      str.alphaX = [str.alphaX, TissueParam{k}.alphaX];
%                  end
% 
%                  if ~isfield(TissueParam{k}, 'betaX') || (isfield(TissueParam{k}, 'beatX') && ~any(ismember(machine.data(1).betaX,TissueParam{k}.betaX)))
%                      %Is matlab evaluating the wholle condition anyway? Does
%                      %it stop if first condition is met?
%                      matRad_cfg.dispWarning(['Warning! Tissue parameter BetaX for tissue ' cst{k,1} 'not available, set do default \n']);
%                      if ismember(obj.default_BetaX, machine.data(1).betaX)
%                          str.betaX = [str.betaX, obj.default_BetaX];
%                      else
%                          str.betaX = [str.betaX, machine.data(1).betaX(1)];
%                          matRad_cfg.dispWarning(['Warning! Tissue parameter default_BetaX for tissue ' cst{k,1} 'not compatible with basedata, setting the available one \n']);
%                      end    
%                  else
%                      str.betaX = [str.betaX, TissueParam{k}.betaX];
%                  end

                 IdxTissue = find(ismember(machine.data(1).alphaX,TissueParam{k}.alphaX) & ...
                           ismember(machine.data(1).betaX,TissueParam{k}.betaX));


                 tissueClass(cst{k,4}{1},1) = IdxTissue;
              end
              str.tissueClass = tissueClass(tissueClass>0);              
           end
       end
       
       function [bixelAlpha,bixelBeta] = calcLQParameter(obj,vRadDepths,baseDataEntry,TissueParam,ix)%mTissueClass,vAlpha_x,vBeta_x,vABratio);
            
            bixelAlpha = NaN*ones(numel(vRadDepths),1);
            bixelBeta  = NaN*ones(numel(vRadDepths),1);
            
            % range shift
            depths = baseDataEntry.depths + baseDataEntry.offset;

            %numOfTissueClass = size(baseDataEntry(1).alpha,2);
             numOfTissueClass = size(unique(TissueParam.tissueClass(ix)),1);
             for i = 1:numOfTissueClass
                 bixelAlpha(TissueParam.tissueClass(ix)==i) = matRad_interp1(depths,baseDataEntry.alpha(:,i),vRadDepths(TissueParam.tissueClass(ix)==i));
                 bixelBeta(TissueParam.tissueClass(ix)==i)  = matRad_interp1(depths,baseDataEntry.beta(:,i), vRadDepths(TissueParam.tissueClass(ix)==i));
             end

       end
   end 
end