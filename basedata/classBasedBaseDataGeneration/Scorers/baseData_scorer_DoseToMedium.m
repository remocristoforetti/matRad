classdef baseData_scorer_DoseToMedium < baseData_genericScorer
    properties
        
    end


    methods (Access = private)
        function obj = baseData_scorer_DoseToMedium()
            obj@baseData_genericScorer();
            obj.name = 'DoseToMedium';
        end
    end

    methods (Static)
        function obj = instance()
            persistent currentInstance;

            if ~isempty(currentInstance)
                obj = currentInstance;
            else
                currentInstance = baseData_scorer_DoseToMedium();
                obj = currentInstance;
            end
        end
    end
end