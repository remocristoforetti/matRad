classdef matRad_bioModel_none < matRad_BiologicalModel
    %% Properities
    properties
        AvailableradiationModalities = {'protons', 'carbon', 'helium'};
        AvailableQuantitiesForOpt = {'physicalDose'};
        RequiredBaseData = {};
    end

    methods
         function obj = matRad_bioModel_none(radiationMode)
         if ~ischar(radiationMode)
           matRad_dispToConsole(['Something wrong with bioModel inputs'],[],'warning');
         end
         
         obj@matRad_BiologicalModel();

         obj.model = 'none';

         obj.radiationMode = radiationMode;
       end
    end
end