classdef matRad_bioModelMKMLET < matRad_bioModelLETBased
%% Properties
    properties (SetAccess = private)
        p0_MKM = 1;
        p1_MKM = 0.229;
        p2_MKM = 1;
        p3_MKM = 0;

        default_RD = 0.72;
        AvailableradiationModalities = {'protons'};
        AvailableQuantitiesForOpt = {'physicalDose','RBExD','effect'};
        RequiredBaseData = {'depths','offset','LET'};
    end

    properties (Hidden)
        
    end
%% Methods
    methods

        function obj = matRad_bioModelMKMLET(radiationMode)
            obj@matRad_bioModelLETBased(radiationMode);
            obj.model = 'MKMLET';
        end

        function p0 = getP0(obj)
            p0 = obj.p0_MKM;
        end

        function p1 = getP1(obj)
            p1 = obj.p1_MKM./((obj.default_RD^2) .* obj.vABratio);
        end

        function p2 = getP2(obj)
            p2 = obj.p2_MKM;
        end

        function p3 = getP3(obj)
            p3 = obj.p3_MKM;
        end
    end
end