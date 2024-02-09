classdef matRad_RobustenssNone < Robustness.matRad_Robustness
    properties
        name = 'none';
    end

    methods
        function this = matRad_RobustenssNone()
            this@Robustness.matRad_Robustness();
        end

        function [d_i,idx]= getDi(this,d,ixContour)

            d_i = cell(numel(this.useScen),1);
            for scenIdx=this.useScen
                idx{scenIdx} = ixContour{this.useScen(scenIdx)};
                d_i{scenIdx} = d{scenIdx}(idx{scenIdx});
            end

        end

        function f_rob = getRobustF(this,f_i)
            f_rob = sum(f_i{:});
        end

        function doseGradient = getRobustGradient(this,doseGradient,doseGrad_i,idxs)


            for scenIdx=1:numel(doseGrad_i)
                currScen = this.useScen(scenIdx);                 
                doseGradient{currScen}(idxs{scenIdx}) = doseGradient{currScen}(idxs{scenIdx}) + doseGrad_i{scenIdx};
            end

        end
        
        function this = assignScenarios(this,multScen)
            % None robustness takes into account all CT scenarios, equally
            % weighted
            nominalCtScenarios = [1:multScen.numOfCtScen];
            this.useScen = nominalCtScenarios;
        end

    end
end