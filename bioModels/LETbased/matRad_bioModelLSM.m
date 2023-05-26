classdef matRad_bioModelLSM < matRad_bioModelLETBased
%% Properties
    properties (SetAccess = private)
        p_lambda_1_1          = 0.008; %0.008; % according to Malte Frese https://www.ncbi.nlm.nih.gov/pubmed/20382482 (FITTED for head and neck patients !)
        p_corrFacEntranceRBE = 0.5;   %[kev/mum]
        p_upperLETThreshold  = 30;    %[kev/mum]
        p_lowerLETThreshold  = 0.3;   %[kev/mum]

        AvailableradiationModalities = {'protons'};
        AvailableQuantitiesForOpt = {'physicalDose','RBExD','effect'};
        RequiredBaseData = {'depths','offset','LET'};
    end

    properties (Hidden)
        
    end
%% Methods
    methods

        function obj = matRad_bioModelLSM(radiationMode)
            obj@matRad_bioModelLETBased(radiationMode);
            obj.model = 'LSM';
        end

        function p0 = getP0(obj)

            %ix = obj.p_lowerLETThreshold < bixelLET < obj.p_upperLETThreshold;
            ix_lower = obj.LET < obj.p_lowerLETThreshold;
            ix_upper = obj.LET > obj.p_upperLETThreshold;

            p0 = obj.vAlpha_x - obj.p_lambda_1_1*obj.p_corrFacEntranceRBE;
            
           
            p0(ix_lower) = p0(ix_lower) + obj.p_lambda_1_1*obj.p_lowerLETThreshold;
            p0(ix_upper) = p0(ix_upper) + obj.p_lambda_1_1*obj.p_upperLETThreshold;
            
            p0 = p0./obj.vAlpha_x;
        end

        function p1 = getP1(obj)
            ix = obj.p_lowerLETThreshold < obj.LET < obj.p_upperLETThreshold;

            p1 = zeros(size(obj.vAlpha_x));
            
            p1(ix) = obj.p_lambda_1_1./obj.vAlpha_x(ix);

        end

        function p2 = getP2(obj)
            p2 = 1;
        end

        function p3 = getP3(obj)
            p3 = 0;
        end
    end
end