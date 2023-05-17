classdef matRad_bioModelCAR < matRad_bioModelLETBased
%% Properties
    properties (SetAccess = private)
        p0_CAR = 0.843;
        p1_CAR = 0.154;
        p2_CAR = 2.686;
        p3_CAR = 1.09;
        p4_CAR = 0.006;

        AvailableradiationModalities = {'protons'};
        AvailableQuantitiesForOpt = {'physicalDose','RBExD','effect'};
        RequiredBaseData = {'depths','offset','LET'};
    end

    properties (Hidden)
        
    end
%% Methods
    methods

        function obj = matRad_bioModelCAR(radiationMode)
            obj@matRad_bioModelLETBased(radiationMode);
            obj.model = 'CAR';
        end

        function p0 = getP0(obj)
            p0 = obj.p0_CAR;
        end

        function p1 = getP1(obj)
            p1 = (obj.p1_CAR * obj.p2_CAR)./obj.vABratio;
        end

        function p2 = getP2(obj)
            p2 = obj.p3_CAR;
        end

        function p3 = getP3(obj)
            p3 = (obj.p4_CAR * obj.p2_CAR)./obj.vABratio;
        end
    end
end