classdef (Abstract) matRad_OptimizationQuantity < handle

    properties (Abstract,Constant)
        quantityName;
        requiredSubquantities;
    end

    properties
        subQuantities;
        d;
        wGrad;
        wCache;
        wGradCache;
    end

    methods
        function this = matRad_OptimizationQuantity()

        end

        function output = getResult(~)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('Function needs to be implemented by subclass');
            output = [];

        end

        function output = getProjectedGradient(~)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('Function needs to be implemented by subclass');
            output = [];
        end

        function output = computeQuantity(~)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('Function needs to be implemented by subclass');
            output = [];
        end

         function output = projectGradient(~)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('Function needs to be implemented by subclass');
            output = [];
        end
    end

    methods (Static)
        function optiFunc = setBiologicalDosePrescriptions(optiFunc,alphaX,betaX)

        end
    end

end