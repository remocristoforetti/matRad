classdef matRad_RobustnessVWWC_INV < Robustness.matRad_Robustness
    properties
        name = 'STOCH';
        %scenProb;
        nominalCTScens;
        
        structureDefinition;
        ix_cache;
    end

    properties (SetAccess= private)

    end
    
    methods
        function this = matRad_RobustnessVWWC_INV()
            this@Robustness.matRad_Robustness();
        end

        function [d_i,idx] = getDi(this,d,ixContour)

            matRad_cfg = MatRad_Config.instance();
            
            idx = ixContour;

            if numel(ixContour)>1
                matRad_cfg.dispError('4D VWWC optimization is currently not supported');
            end
            
            if ~exist('d_tmp','var')
                d_tmp = [d{this.useScen}];
            end

            d_Scen = d_tmp(ixContour{1},:);

            if isequal(this.structureDefinition,'OAR')
                [d_min, ix] = min(d_Scen,[],2);
                d_i = {d_min};
            elseif isequal(this.structureDefinition,'TARGET')
                [d_max, ix] = max(d_Scen,[],2);
                d_i = {d_max};
            end

            if ~isempty(this.ix_cache)
                matRad_cfg.dispError('This should not happen');
            else
                this.ix_cache = ix;
            end
        end

        function f_rob = getRobustF(this,f_i)
            f_rob = f_i{1};
            this.ix_cache = [];
        end

        function doseGradient = getRobustGradient(this,doseGradient,doseGrad_i,idx)

            matRad_cfg = MatRad_Config.instance();
            
             if isempty(this.ix_cache)
                matRad_cfg.dispError('This should not happen');
            else
                ix = this.ix_cache;
             end
            

             for scenIdx = 1:numel(this.useScen)
  
                currScen = this.useScen(scenIdx);

                currWcIx = double(ix == scenIdx);
    
                % 4D not supported for VWWC
                doseGradient{currScen}(idx{1}) = doseGradient{currScen}(idx{1}) + doseGrad_i{1}.*currWcIx;
            end

            this.ix_cache = [];

        end
        
        function this = assignScenarios(this,multScen)
            this.useScen        = find(multScen.scenMask);
            %keep them all for the time being
            this.nominalCTScens = multScen.linearMask(:,1);
        end

        function set.structureDefinition(this,value)
            
            matRad_cfg = MatRad_Config.instance();
            valid = ischar(value) && any(strcmp(value, {'OAR', 'TARGET'}));

            if valid
                this.structureDefinition = value;
            else
                matRad_cfg.dispError('Unable to set structure defintion %', value);
            end

        end
    end
end