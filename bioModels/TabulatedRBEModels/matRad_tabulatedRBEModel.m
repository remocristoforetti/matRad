classdef matRad_tabulatedRBEModel < matRad_BiologicalModel
   %bioOpt is allowed even if there is no table, should correct for that.
   %-> maybe in constructor turn it off
   properties
      %Do these make sense here?

      default_alphaX = 0.1;
      default_betaX = 0.05;

      RBEtable;
      weightMode = 'eD';
   end
   
   properties
       AvailableradiationModalities = {'protons', 'carbon'};
       AvailableQuantitiesForOpt = {'physicalDose','RBExD','effect'};
       RequiredBaseData = {'depths','offset', 'RBEtable'};
      
   end

   methods

      function obj = matRad_tabulatedRBEModel(radiationMode)
         
         obj@matRad_BiologicalModel();
         obj.radiationMode = radiationMode;
         %obj.bioOpt = false;
      end
      
      
      function set.RBEtable(obj,RBEtableName)
         
         load([RBEtableName, '.mat']);
         obj.model = RBEtable.meta.model;
         obj.RBEtable = RBEtableName;
         obj.bioOpt = true;
      end

      function str = calcTissueParameters(obj,cst,numVoxels,stf,ctScen)
         %alpha/beta table should have:
         %table.meta.model
         
         %table.data(n).energy; -> cell array for different ions
         
         %table.data(n).alpha; -> matrix ExIons
         %table.data(n).beta;
         %table.data(n).alphaX;
         %table.data(n).betaX;
         %and different alphaX/betaX go into different n
         
         %What is generic for this class is calcLQParameters -> given
         %varAlpha/varBeta as function of depth, just sample it.
         
         %calcTissueParameters is specific of submodule
         matRad_cfg = MatRad_Config.instance();
         str = struct('alphaX', [], ...
                      'betaX', [], ...
                      'tissueABratio', [], ...
                      'tissueClass', [], ...
                      'baseParam', []);
         tissueParam = {cst{:,5}};
          
         tissueClass = zeros(numVoxels,1);
         for tisueIdx = 1:size(tissueParam,2)
              
            if ~isfield(tissueParam{tisueIdx}, 'alphaX')
                str.alphaX = [str.alphaX, obj.default_AlphaX];
                matRad_cfg.dispWarning(['Warning! Tissue parameter AlphaX for tissue ' cst{tisueIdx,1} 'not available, set do default \n']);
            else
                str.alphaX = [str.alphaX, tissueParam{tisueIdx}.alphaX];
            end
              
            if ~isfield(tissueParam{tisueIdx}, 'betaX')
                str.betaX = [str.betaX, obj.default_BetaX];
                matRad_cfg.dispWarning(['Warning! Tissue parameter BetaX for tissue ' cst{tisueIdx,1} 'not available, set do default \n']);
            else
                str.betaX = [str.betaX, tissueParam{tisueIdx}.betaX];
            end
              
            str.tissueABratio(tisueIdx) = str.alphaX(tisueIdx)./str.betaX(tisueIdx);
            tissueClass(cst{tisueIdx,4}{1},1) = tisueIdx; %What if there is an overlap between structures?
         end
         
         str.tissueClass = tissueClass(tissueClass>0);
         
         load([obj.RBEtable, '.mat']);
         
         
         tabulatedAlphaX = [RBEtable.meta.modelParameters.alphaX];
         tabulatedBetaX  = [RBEtable.meta.modelParameters.betaX];
         
         for tissueIdx = 1:size(tissueParam,2)
            alphaLogical = str.alphaX(tissueIdx) == tabulatedAlphaX;
            betaLogical = str.betaX(tissueIdx) == tabulatedBetaX;
            RBEtableIdx = find(alphaLogical & betaLogical);
            RBEdata(tissueIdx) = RBEtable.data(RBEtableIdx);

         end

         energiesInitialBeam = unique([stf.ray.energy]);
         %load energies from DS energy binning
         try 
            load([obj.radiationMode, '_', obj.baseData]);
         catch
            matRad_cfg.dispError('Unable to load machine data %s',[obj.radiationMode, '_', obj.machine]);
         end
         
         [~,eIdx] = intersect([machine.data(:).energy],energiesInitialBeam);

         switch obj.radiationMode
            case 'protons'
               nFragments = 1;
            case 'carbon'
               nFragments = 6;
         end
         

         for beamEnergyIdx=1:size(eIdx,1) 
         
             tmp1 = [];

            switch obj.weightMode
                  case 'eD'
                     Phi = machine.data(eIdx(beamEnergyIdx)).DS.edPhi;%eDPhi;
                     E  = machine.data(eIdx(beamEnergyIdx)).DS.eBinning;
                     for ionIdx=1:nFragments
                        tmp1(:,:,ionIdx) = full(Phi{ionIdx})';
                     end
                  case 'LET'
                     Phi = machine.data(eIdx(beamEnergyIdx)).DS.Phi;
                     E  = machine.data(eIdx(beamEnergyIdx)).DS.eBinning;
                     for ionIdx=1:nFragments
                        let(:,ionIdx) = obj.computeLET(E,ionIdx)';
                        tmp1(:,:,ionIdx) = (full(Phi{ionIdx}).*repmat(let(:,ionIdx),1,size(Phi{ionIdx},2)))';
                     end
             end
            
            %Phi = machine.data(eIdx(beamEnergyIdx)).DS.Phi;
            %E  = machine.data(eIdx(beamEnergyIdx)).DS.EBinning;
            %Ed = machine.data(eIdx(beamEnergyIdx)).DS.EdBinning;
            
            

            %varNumAlpha = 0;
            %varNumBeta = 0;
            %varDen = 0;

            str.baseParam(beamEnergyIdx).varAlpha = zeros(size(Phi{1},2),size(tissueParam,2));
            str.baseParam(beamEnergyIdx).varBeta = zeros(size(Phi{1},2),size(tissueParam,2));
            
            for tissueIdx=1:size(tissueParam,2)
               currRBEdata = RBEdata(tissueIdx);

               alphaE = interp1(currRBEdata.energies,currRBEdata.alpha(:,1:nFragments),E');
               betaE = interp1(currRBEdata.energies,currRBEdata.beta(:,1:nFragments),E');
               
               alphaE = permute(alphaE, [1,3,2]);
               betaE = permute(betaE, [1,3,2]);
               
              
               %Here sum is carried out on number of fragments;
               varAlphaNum = squeeze(sum(pagemtimes(tmp1,alphaE), [3]));
               varDen = squeeze(sum(tmp1,[2,3]));               
               varBetaNum = squeeze(sum(pagemtimes(tmp1,betaE), [3]));

               Idx = varDen > 0;
               
               str.baseParam(beamEnergyIdx).beamEnergy = machine.data(eIdx(beamEnergyIdx)).energy;
               str.baseParam(beamEnergyIdx).varAlpha(Idx,tissueIdx) = varAlphaNum(Idx)./varDen(Idx);
               str.baseParam(beamEnergyIdx).varBeta(Idx,tissueIdx) = varBetaNum(Idx)./varDen(Idx);
               
               %matRad_cfg.dispWarning('!!!! 3 mm shift inserted in tabulatedModel !!!');
               str.baseParam(beamEnergyIdx).d = machine.data(eIdx(beamEnergyIdx)).DS.depths;% + 3;
               %end
            
            end
            %alphaE = interp1(RBEtable.)
                       
         end
      end
      
      function [bixelAlpha, bixelBeta] = calcLQParameter(obj,vRadDepths,baseDataEntry,TissueParam,ix)
          bixelAlpha = NaN*ones(numel(vRadDepths),1);
          bixelBeta  = NaN*ones(numel(vRadDepths),1);
          
          % For now
          [~,eIdx] = intersect([TissueParam.baseParam.beamEnergy], baseDataEntry.energy);
          depths = TissueParam.baseParam(eIdx).d;
          
          numOfTissueClass = size(unique(TissueParam.tissueClass),1);
          for i=1:numOfTissueClass
             currIdx = TissueParam.tissueClass(ix)==i;
             bixelAlpha(currIdx) = interp1(depths,TissueParam.baseParam(eIdx).varAlpha(:,i), vRadDepths(currIdx));
             bixelBeta(currIdx) = interp1(depths,TissueParam.baseParam(eIdx).varBeta(:,i), vRadDepths(currIdx));
          end
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
   end
   
end