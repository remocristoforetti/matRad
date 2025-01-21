classdef matRad_PhysicalDoseProtons < matRad_PhysicalDoseSingleModality

    properties (Constant)
        quantityName = 'physicalDoseProtons';
        requiredSubquantities = {};
        
        modality = 'protons';
    end

    methods
        
        function this = matRad_PhysicalDoseProtons()
            this@matRad_PhysicalDoseSingleModality;
        end
    end
end