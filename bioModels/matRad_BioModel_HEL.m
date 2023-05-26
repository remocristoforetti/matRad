classdef matRad_bioModel_HEL < matRad_Bio_Model
      properties (SetAccess = private)%constant
        
        p1_HEL = NaN;
        p2_HEL = NaN;
        p0_HEL = NaN;
        default_alphaX = 0.1;
        default_betaX = 0.05;

        AvailableradiationModalities = {'helium'};
        AvailableQuantitiesForOpt = {'physicalDose','RBExD','effect'};
        RequiredBaseData = {'depths','offset','LET'};
        NumModelParameters = 2;
        ModelParameters = {'alphaX','betaX'};
      end


      %% Methods
   methods
       %Constructor
       function obj = matRad_bioModel_HEL(radiationMode)
          if ~ischar(radiationMode)
           matRad_dispToConsole(['Something wrong with bioModel inputs'],[],'warning');
          end
          
          obj@matRad_Bio_Model();
%          obj.model = obj.Model_name;
          obj.model = 'HEL';
          obj.radiationMode = radiationMode;
       end
       
       function str = calcTissueParameters(obj,cst,numVoxels,ctScen)
           matRad_cfg =  MatRad_Config.instance();

           str = struct('alphaX', [], ...
                        'betaX', [], ...
                        'tissueABratio', [], ...
                        'tissueClass', []);
           
           TissueParam = {cst{:,5}};
           tissueClass = zeros(numVoxels,1);
           
           try
              load([matRad_cfg.matRadRoot filesep 'basedata' filesep strcat(obj.radiationMode, '_', obj.baseData) '.mat']);
           catch
              matRad_cfg.dispError(['Could not find the following machine file: ' obj.machine '\n']);
           end

           if any(machine)
              for k = 1:size(TissueParam,2)

                 if ~isfield(TissueParam{k}, 'alphaX') || (isfield(TissueParam{k}, 'alphaX') && ~any(ismember(machine.data(1).alphaX,TissueParam{k}.alphaX)))
                     matRad_cfg.dispWarning(['Warning! Tissue parameter AlphaX for tissue ' cst{k,1} 'not available, set do default \n']);
                     if ismember(obj.default_AlphaX, machine.data(1).alphaX)
                         str.alphaX = [str.alphaX, obj.default_AlphaX];
                     else
                         str.alphaX = [str.alphaX, machine.data(1).alphaX(1)];
                         matRad_cfg.dispWarning(['Warning! Tissue parameter default_AlphaX for tissue ' cst{k,1} 'not compatible with basedata, setting the available one \n']);
                     end
                 else
                     str.alphaX = [str.alphaX, TissueParam{k}.alphaX];
                 end

                 if ~isfield(TissueParam{k}, 'betaX') || (isfield(TissueParam{k}, 'beatX') && ~any(ismember(machine.data(1).betaX,TissueParam{k}.betaX)))
                     matRad_cfg.dispWarning(['Warning! Tissue parameter BetaX for tissue ' cst{k,1} 'not available, set do default \n']);
                     if ismember(obj.default_BetaX, machine.data(1).betaX)
                         str.betaX = [str.betaX, obj.default_BetaX];
                     else
                         str.betaX = [str.betaX, machine.data(1).betaX(1)];
                         matRad_cfg.dispWarning(['Warning! Tissue parameter default_BetaX for tissue ' cst{k,1} 'not compatible with basedata, setting the available one \n']);
                     end    
                 else
                     str.betaX = [str.betaX, TissueParam{k}.betaX];
                 end
                 
                 str.tissueABratio(k) = str.alphaX(k)./str.betaX(k);
                 tissueClass(cst{k,4}{1},1) = k;
              end
              str.tissueClass = tissueClass(tissueClass>0);
           end
       end

       function [bixelAlpha,bixelBeta] = calcLQParameter(obj,vRadDepths,baseDataEntry,TissueParam,ix);
            
            bixelAlpha = NaN*ones(numel(vRadDepths),1);
            bixelBeta  = NaN*ones(numel(vRadDepths),1);
            
            % range shift
            depths = baseDataEntry.depths + baseDataEntry.offset;
 
                              % data-driven RBE parametrization of helium ions
                  % https://iopscience.iop.org/article/10.1088/0031-9155/61/2/888
                  
                  bixelLET = matRad_interp1(depths,baseDataEntry.LET,vRadDepths);
                  bixelLET(isnan(bixelLET)) = 0;
                  
                  % quadratic fit
                  %f_Q      = 8.53959e-4 .* bixelLET.^2;
                  %RBEmax_Q = 1 + 2.145e-1  + vABratio.^-1 .* f_Q;
                  % linear quadratic fit
                  %f_LQ      = 2.91783e-1*bixelLET - 9.525e-4*bixelLET.^2;
                  %RBEmax_LQ = 1 + ((1.42057e-1 + (vABratio.^-1)) .* f_LQ);
                  % linear exponential fit
                  %f_LE      = (2.965e-1 * bixelLET) .* exp(-4.90821e-3 * bixelLET);
                  %RBEmax_LE = 1 + ((1.5384e-1  + (vABratio.^-1)) .* f_LE);
                  vAlpha_x = TissueParam.alphaX(TissueParameter.tissueClass(ix));
                  vBeta_x = TissueParam.betaX(TissueParameter.tissueClass(ix));
                  vABRatio = TissueParam.tissueABratio(TisseuParameter.tissueClass(ix));
                  
                  % quadratic exponential fit
                  f_QE      = (obj.p1_HEL * bixelLET.^2) .* exp(-obj.p2_HEL * bixelLET);
                  RBEmax_QE = 1 + ((obj.p0_HEL  + (vABratio.^-1)) .* f_QE);

                  % the linear quadratic fit yielded the best fitting result
                  RBEmax = RBEmax_QE;
                  
                  RBEmin = 1; % no gain in using fitted parameters over a constant value of 1
                  
                  bixelAlpha = RBEmax    .* vAlpha_x;
                  bixelBeta  = RBEmin.^2 .* vBeta_x;
                    
       end
   end 
end