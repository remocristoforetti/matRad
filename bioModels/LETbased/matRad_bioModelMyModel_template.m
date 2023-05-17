classdef matRad_bioModelMyModel_template < matRad_bioModelLETBased
% LET based Biological model for protons.
% The model computes the variable alpha and beta parameters for protons
% following:
    %  bixelAlpha = RBEmax.*AlphaX;
    %  bixelBeta = (RBEmin.^2).*BetaX;

% where:
    %  RBEmax = p0 + p1.*LET;
    %  RBEmin = p2 + p3.*LET;
% and AlphaX and BetaX refer to the alpha/beta parameters for the reference
% radiation



%% Properties
    properties (SetAccess = private)

        %Define the model parameters
        p0_MyModel_template = 1;         
        p1_MyModel_template = 1;
        p2_MyModel_template = 1;
        p3_MyModel_template = 1;

        AvailableradiationModalities = {'protons'};
        AvailableQuantitiesForOpt = {'physicalDose','RBExD','effect'};
        RequiredBaseData = {'depths','offset','LET'};
    end

    properties (Hidden)
        
    end
%% Methods
    methods

        function obj = matRad_bioModelMyModel_template(radiationMode)
            obj@matRad_bioModelLETBased(radiationMode);
            obj.model = 'MyModel';
        end

        % These methods are provided for specific LET-based model
        % implementation
        function p0 = getP0(obj)
            p0 = obj.p0_MyModel_template;
        end

        function p1 = getP1(obj)
            p1 = obj.p1_MyModel_template./obj.vABratio;
        end

        function p2 = getP2(obj)
            p2 = obj.p2_MyModel_template;
        end

        function p3 = getP3(obj)
            p3 = obj.p3_MyModel_template;
        end
    end
end