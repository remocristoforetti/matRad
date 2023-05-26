classdef matRad_bioModel_LEM < matRad_BiologicalModel

    properties (SetAccess = private)
        AvailableradiationModalities = {'carbon', 'helium'};                        %radiation modalities that can be optimized with this model
        AvailableQuantitiesForOpt = {'physicalDose','RBExD','effect'};              %optimization quantities that can be optimized with this model
        RequiredBaseData = {'depths','offset','alpha', 'beta', 'alphaX', 'betaX'};  %required fields in basData
        default_AlphaX = 0.1;
        default_BetaX = 0.05;
   end

      %% Methods

      methods
       %Constructor
       function obj = matRad_bioModel_LEM(radiationMode)
          if ~ischar(radiationMode)
           matRad_cfg.dispWarning(['Something wrong with bioModel inputs']);
          end
          
          obj@matRad_BiologicalModel();
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
            depths = baseDataEntry.depths;% + baseDataEntry.offset;

            %numOfTissueClass = size(baseDataEntry(1).alpha,2);
             numOfTissueClass = size(unique(TissueParam.tissueClass),1);
             tissueIdxs = unique(TissueParam.tissueClass);
             for i = 1:numOfTissueClass
                 bixelAlpha(TissueParam.tissueClass(ix)==tissueIdxs(i)) = matRad_interp1(depths,baseDataEntry.alpha(:,i),vRadDepths(TissueParam.tissueClass(ix)==tissueIdxs(i)));
                 bixelBeta(TissueParam.tissueClass(ix)==tissueIdxs(i))  = matRad_interp1(depths,baseDataEntry.beta(:,i), vRadDepths(TissueParam.tissueClass(ix)==tissueIdxs(i)));
             end

       end
   end 
end