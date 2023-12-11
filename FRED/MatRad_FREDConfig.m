classdef MatRad_FREDConfig < handle

    properties
        patientFilename      = 'CTpatient.mhd';
        runInputFilename     = 'fred.inp';
        
        regionsFilename      = 'regions.inp';
        funcsFilename        = 'funcs.inp';
        planFilename         = 'plan.inp';
        fieldsFilename       = 'fields.inp';
        layersFilename       = 'layers.inp';
        beamletsFilename     = 'beamlets.inp';
        planDeliveryFilename = 'planDelivery.inp';

        hLutLimits = [-1000,1375];
        conversionFactor = 1e6;
        planDeliveryTemplate = 'planDelivery.txt';
        FREDrootFolder;

        MCrunFolder;
        inputFolder;
        regionsFolder;
        planFolder;

    end

    methods (Access = private)
        function obj = MatRad_FREDConfig()
        end
    end

    methods (Static)
        function obj = instance()
            %instance creates a singleton instance
            %  In MatRad_FREDConfig, the constructor is private to make sure only on global instance exists.
            %  Call this static function to get or create an instance of the matRad_FRED configuration class
            persistent uniqueInstance;
            
            if isempty(uniqueInstance)
                obj = MatRad_FREDConfig();
                uniqueInstance = obj;
            else
                obj = uniqueInstance;
            end
        end
    
    end

    methods
        function set.FREDrootFolder(obj, pathValue)
            obj.FREDrootFolder = pathValue;

            obj.updatePaths;
        end


        function updatePaths(obj)
            obj.MCrunFolder     = fullfile(obj.FREDrootFolder, 'MCrun');
            obj.inputFolder     = fullfile(obj.MCrunFolder, 'inp');
            obj.regionsFolder   = fullfile(obj.inputFolder, 'regions');
            obj.planFolder      = fullfile(obj.inputFolder, 'plan');
        end
    end



end