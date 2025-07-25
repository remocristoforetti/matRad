classdef (Abstract) matRad_LexicographicIteration < handle

    properties
        slack;
    end

    properties (Constant)
        defaultSlack =  1.2;
    end

    methods

        function this = matRad_LexicographicIteration(pln)

            if nargin>0
                if isfield(pln, 'propOpt') && isfield(pln.propOpt, 'slack')
                    this.slack = pln.propOpt.slack;
                else
                    this.slack = this.defaultSlack;
                end
            else
                this.slack = this.defaultSlack;
            end
        end

        function [dij,cst,pln,wInit,optiProb, meta] = getProblem(this,list, dij, originalCst, pln, wInit)
            
            % scroll the list and get the next objecive to be optimized
            [objective, currPriority] = this.scrollList(list);

            % build the cst form list for optimization
            cst = this.getCst(list, originalCst, objective, currPriority);

            % Init the optimization/quantities and so on
            [dij,cst,pln,wInit,optiProb] = matRad_initOptimization(dij,cst,pln,wInit);

            % Additional meta information to be passed on to the lexEngine
            meta = this.addProblemInfo(objective, currPriority,dij,cst,pln,wInit,optiProb);
            
            meta.stepPriority = currPriority;
        end


        function [objective, priority] = scrollList(this,~)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('This function needs to be implemented by the subclass');
        end

        function cst = getCst(this,~, ~, ~, ~)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('This function needs to be implemented by the subclass');
        end

        function meta = addProblemInfo(~,~,~,~,~,~,~,~)
            % Additional function that is called at the end of tge
            % getProblem method and allows custom information to be passed
            % on in the meta structure to the lexicographic engine.
            meta = [];
        end
    end
end