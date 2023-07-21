classdef matRad_bioModelMCN < matRad_bioModelLETBased
% subclass that implements the MCN model
% (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4634882/) (accessed on 21/7/2023)
%% Properties
    properties (SetAccess = private)
        p0_MCN = 0.999064;
        p1_MCN = 0.35605;
        p2_MCN = 1.1012;
        p3_MCN = -0.0038703;

        AvailableRadiationModalities = {'protons'};
    end

    properties (Hidden)
        
    end
%% Methods
    methods

        function obj = matRad_bioModelMCN(radiationMode)
            obj@matRad_bioModelLETBased(radiationMode);
            obj.model = 'MCN';
        end

        function p0 = getP0(obj)
            p0 = obj.p0_MCN;
        end

        function p1 = getP1(obj)
            p1 = obj.p1_MCN./obj.vABratio;
        end

        function p2 = getP2(obj)
            p2 = obj.p2_MCN;
        end

        function p3 = getP3(obj)
            p3 = obj.p3_MCN.*sqrt(obj.vABratio);
        end
    end
end