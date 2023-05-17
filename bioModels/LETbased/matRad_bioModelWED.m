classdef matRad_bioModelWED < matRad_bioModelLETBased
%% Properties
    properties (SetAccess = private)
        p0_WED = 1;
        p1_WED = 0.434;
        p2_WED = 1;
        p3_WED = 0;

        AvailableradiationModalities = {'protons'};
        AvailableQuantitiesForOpt = {'physicalDose','RBExD','effect'};
        RequiredBaseData = {'depths','offset','LET'};
    end

    properties (Hidden)
        
    end
%% Methods
    methods

        function obj = matRad_bioModelWED(radiationMode)
            obj@matRad_bioModelLETBased(radiationMode);
            obj.model = 'WED';
        end

        function p0 = getP0(obj)
            p0 = obj.p0_WED;
        end

        function p1 = getP1(obj)
            p1 = obj.p1_WED./obj.vABratio;
        end

        function p2 = getP2(obj)
            p2 = obj.p2_WED;
        end

        function p3 = getP3(obj)
            p3 = obj.p3_WED;
        end
    end
end