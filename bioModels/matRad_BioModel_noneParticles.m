classdef matRad_BioModel_noneParticles < matRad_Bio_Model
    %% Properties
    properties (SetAccess = private)%constant

        AvailableradiationModalities = {'protons','carbon','helium'};
        AvailableQuantitiesForOpt = {'physicalDose'};
        RequiredBaseData = {};

    end 
    %% Methods
   methods
       function obj = matRad_BioModel_noneParticles(radiationMode)
            if ~ischar(radiationMode)
               matRad_cfg.dispWarning(['Something wrong with bioModel inputs']);
            end
            
           
            obj@matRad_Bio_Model();
            obj.model = 'none';
            obj.radiationMode = radiationMode;
       end
   end
end