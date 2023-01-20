classdef matRad_BioModel_LSM < matRad_Bio_Model
    
       properties (SetAccess = private)%constant
        p_lamda_1_1          = 0.008; %0.008; % according to Malte Frese https://www.ncbi.nlm.nih.gov/pubmed/20382482 (FITTED for head and neck patients !)
        p_corrFacEntranceRBE = 0.5;   %[kev/mum]
        p_upperLETThreshold  = 30;    %[kev/mum]
        p_lowerLETThreshold  = 0.3;   %[kev/mum]
        
        AvailableradiationModalities = {'protons'};
        AvailableQuantitiesForOpt = {'physicalDose','RBExD','effect'};
        RequiredBaseData = {'depths','offset','LET'};
        NumModelParameters = 2;
        ModelParameters = {'alphaX','betaX'};
   end

      %% Methods
   methods
       %Constructor
       function obj = matRad_BioModel_LSM(radiationMode)
          if ~ischar(radiationMode)
           matRad_cfg.dispWarning(['Something wrong with bioModel inputs']);
          end
          
          obj@matRad_Bio_Model();
          %obj.model = obj.Model_name;
          obj.model = 'LSM';
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
                  matRad_dispToConsole(['Warning! Tissue parameter AlphaX for tissue ' cst{k,1} 'not available, set do default \n'],[],'warning');
              else
                  str.alphaX = [str.alphaX, TissueParam{k}.alphaX];
              end
               
              if ~isfield(TissueParam{k}, 'betaX')
                  str.betaX = [str.betaX, obj.default_BetaX];
                  matRad_dispToConsole(['Warning! Tissue parameter BetaX for tissue ' cst{k,1} 'not available, set do default \n'],[],'warning');
              else
                  str.betaX = [str.betaX, TissueParam{k}.betaX];
              end
              
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

            ixd                 = obj.p_lowerLETThreshold < bixelLET < obj.p_upperLETThreshold;

            alpha_0            = vAlpha_x - (obj.p_lamda_1_1 * obj.p_corrFacEntranceRBE);
            bixelAlpha(ixd)  = alpha_0(ixd) + obj.p_lamda_1_1 * bixelLET;

            if sum(ixd) < length(bixelLET)
                bixelAlpha(bixelLET > pln.bioParam.lowerLETThreshold) =  alpha_0(bixelLET > obj.p_upperLETThreshold) + obj.p_lamda_1_1 * obj.p_upperLETThreshold;
                bixelAlpha(bixelLET < pln.bioParam.lowerLETThreshold) =  alpha_0(bixelLET < obj.p_lowerLETThreshold) + obj.p_lamda_1_1 * obj.p_lowerLETThreshold;
            end
            bixelBeta        = vBeta_x;
       end
   end
    
    

end