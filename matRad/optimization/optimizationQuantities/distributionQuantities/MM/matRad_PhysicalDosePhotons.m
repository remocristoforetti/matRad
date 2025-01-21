classdef matRad_PhysicalDosePhotons < matRad_PhysicalDoseSingleModality

    properties (Constant)
        quantityName = 'physicalDosePhotons';
        requiredSubquantities = {};
        
        modality = 'photons';
    end

    methods
        
        function this = matRad_PhysicalDosePhotons()
            this@matRad_PhysicalDoseSingleModality;
        end
    end
end