classdef matRad_bioModelLETBased < matRad_BiologicalModel

    properties
        defaultAlphaX = 0.1;
        defaultBetaX  = 0.05;
    end

    properties (Hidden)
        vABratio;
        p0;
        p1;
        p2;
        p3;
    end
    
    methods
        function obj = matRad_bioModelLETBased(radiationMode)
            obj@matRad_BiologicalModel();
            obj.radiationMode = radiationMode;
        end

        function str = calcTissueParameters(obj,cst,numVoxels,stf,ctScen)
            str = struct('alphaX', [], ...
                        'betaX', [], ...
                        'tissueABratio', [], ...
                        'tissueClass', []);
           
           tissueParam = {cst{:,5}};
           tissueClass = zeros(numVoxels,1);
           for k = 1:size(tissueParam,2)
               
              if ~isfield(tissueParam{k}, 'alphaX')
                  str.alphaX = [str.alphaX, obj.defaultAlphaX];
                  matRad_cfg.dispWarning(['Warning! Tissue parameter AlphaX for tissue ' cst{k,1} 'not available, set do default \n']);
              else
                  str.alphaX = [str.alphaX, tissueParam{k}.alphaX];
              end
               
              if ~isfield(tissueParam{k}, 'betaX')
                  str.betaX = [str.betaX, obj.defaultBetaX];
                  matRad_cfg.dispWarning(['Warning! Tissue parameter BetaX for tissue ' cst{k,1} 'not available, set do default \n']);
              else
                  str.betaX = [str.betaX, tissueParam{k}.betaX];
              end
              
              str.tissueABratio(k) = str.alphaX(k)./str.betaX(k);
              tissueClass(cst{k,4}{1},1) = k; %What if there is an overlap between structures?
           end
           str.tissueClass = tissueClass(tissueClass>0);

        end

        function [bixelAlpha,bixelBeta] = calcLQParameter(obj,vRadDepths,baseDataEntry,tissueParam,ix)
            
            bixelAlpha = NaN*ones(numel(vRadDepths),1);
            bixelBeta  = NaN*ones(numel(vRadDepths),1);

            depths = baseDataEntry.depths + baseDataEntry.offset;

             bixelLET = matRad_interp1(depths,baseDataEntry.LET,vRadDepths);
             bixelLET(isnan(bixelLET)) = 0;

             vAlpha_x = tissueParam.alphaX(tissueParam.tissueClass(ix))';
             vBeta_x = tissueParam.betaX(tissueParam.tissueClass(ix))';
             obj.vABratio = tissueParam.tissueABratio(tissueParam.tissueClass(ix))';

             RBEmax = obj.p0 + obj.p1.*bixelLET;
             RBEmin = obj.p2 + obj.p3.*bixelLET;

             bixelAlpha = RBEmax.*vAlpha_x;
             bixelBeta = (RBEmin.^2).*vBeta_x;
        end

        %Shadow the get functions so that they could be implemented by the
        %subclass
        function p0 = get.p0(obj)
            p0 = obj.getP0();
        end

        function p1 = get.p1(obj)
            p1 = obj.getP1();
        end

        function p2 = get.p2(obj)
            p2 = obj.getP2();
        end

        function p3 = get.p3(obj)
            p3 = obj.getP3();
        end

    end

    methods (Abstract)
       getP0();

       getP1();

       getP2();

       getP3();
    end
end