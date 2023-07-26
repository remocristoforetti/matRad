classdef baseData_phantom_cylinder < baseData_genericPhantom
    properties
        HL;
        rMax;
        rMin;
        rBins;
        zBins;
        material;
        parent;
    end


    methods

        function obj = baseData_phantom_cylinder(name,HL,rMax,rMin,material, varargin)
            
            pars = inputParser();

            addOptional(pars, 'rBins', 1);
            addOptional(pars, 'zBins', 1);
            addOptional(pars, 'parent', 'World');
            parse(pars, varargin{:});

            obj@baseData_genericPhantom();

            obj.name     = name;
            obj.HL       = HL;
            obj.rMax     = rMax;
            obj.rMin     = rMin;
            obj.rBins    = pars.Results.rBins;
            obj.zBins    = pars.Results.zBins;
            obj.parent   = pars.Results.parent;
            obj.material = material;
        end


        function writePhantomParameters(obj,fID)
            %Write code to output the topas lines on fID. This should be
            %called by phantom_setup in a for loop over the phantoms
        end

        function writeScorers(obj,fID)
            %This is called by the scoredQuantity, on the same scorer file
            %need to have all the scorers (from different phantoms) that
            %score same quantity
            if ~isempty(obj.scorers)
                for scorerIdx=1:length(obj.scorers)
                    obj.scorer(scorerIdx).writeScorerParameters(fID,obj.name);
                end
            end
        end
    end
end