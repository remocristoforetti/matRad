classdef matRad_RobustnessOWC < Robustness.matRad_Robustness
    properties
        name = 'STOCH';
        %scenProb;
        nominalCTScens;
        
    end

    properties (SetAccess= private)
        fGrad;
    end
    
    methods
        function this = matRad_RobustnessOWC()
            this@Robustness.matRad_Robustness();
        end

        function [d_i,idx] = getDi(this,d,ixContour)

            d_i = cell(numel(this.useScen),1);
            for scenIdx=1:numel(this.useScen)
                idx{scenIdx} = ixContour{this.nominalCTScens(scenIdx)};
                d_i{scenIdx} = d{this.useScen(scenIdx)}(idx{scenIdx});
            end
        end

        function f_rob = getRobustF(this,f_i)
            f_rob = f_i;
        end


        function doseGradient = getRobustGradient(this,doseGradient,doseGrad_i,idx)

            matRad_cfg = MatRad_Config.instance();

            for scenIdx = 1:numel(this.useScen)
                currScen = this.useScen(scenIdx);

                if this.fGrad(scenIdx) ~= 0
                    doseGradient{currScen}(idx{scenIdx}) = doseGradient{currScen}(idx{scenIdx}) + this.fGrad(scenIdx)*doseGrad_i{scenIdx};
                end
            end

        end
        
        function this = assignScenarios(this,multScen)
            this.useScen        = find(multScen.scenMask);
            %keep them all for the time being
            this.nominalCTScens = multScen.linearMask(:,1);
        end

        function updatefGrad(this, value)
            this.fGrad = value;
        end
    end
end