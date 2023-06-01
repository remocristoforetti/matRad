classdef matRad_bioModel_none < matRad_BiologicalModel
    
    properties
         AvailableRadiationModalities = {'photons','protons', 'carbon', 'helium'};
         RequiredBaseData = {};
    end

    methods
        function obj = matRad_bioModel_none(radiationMode)
            obj@matRad_BiologicalModel(radiationMode);
            obj.model = 'none';
        end
    end
end