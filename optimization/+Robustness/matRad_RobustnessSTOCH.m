classdef matRad_RobustnessSTOCH < Robustness.matRad_Robustness
    properties
        name = 'STOCH';
        scenProb;
        nominalCTScens;
    end

    
    methods
        function this = matRad_RobustnessSTOCH()
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
            f_rob = [f_i{:}]*this.scenProb;
        end

        function doseGradient = getRobustGradient(this,doseGradient,doseGrad_i,idxs)

            doseGrad_i = arrayfun(@(idx)doseGrad_i{idx}*this.scenProb(idx), [1:numel(this.useScen)]', 'UniformOutput', false);

            for scenIdx=1:numel(doseGrad_i)
                currScen = this.useScen(scenIdx);
                doseGradient{currScen}(idxs{scenIdx}) = doseGradient{currScen}(idxs{scenIdx}) + doseGrad_i{scenIdx};
            end

        end
        
        function this = assignScenarios(this,multScen)
            this.useScen        = find(multScen.scenMask);
            this.scenProb       = multScen.scenWeight;
            %keep them all for the time being
            this.nominalCTScens = multScen.linearMask(:,1);
        end
    end
end