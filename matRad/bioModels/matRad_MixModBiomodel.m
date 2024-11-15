classdef matRad_MixModBiomodel < matRad_BiologicalModel

    properties (Constant)
        model = 'MixMod'
        possibleRadiationModes = {'MixMod'};
        requiredQuantities = {};
        defaultReportQuantity = 'physicalDose';
    end

    properties
        singleModalityModels;
    end

    methods
        function this = matRad_MixModBiomodel()
            this = this@matRad_BiologicalModel();
        end
    end
end