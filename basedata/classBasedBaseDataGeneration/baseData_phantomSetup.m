classdef baseData_phantomSetup < handle

    properties
        nPhantoms;
    end

    properties (SetAccess = protected)
        phantoms;
    end

    methods
        function obj = baseData_phantomSetup()
            obj.nPhantoms = 0;
        end

        function addPhantom(obj, phantom)
            matRad_cfg = MatRad_Config.instance();

            currentPhantomNames = arrayfun(@(phantomIdx) obj.phantoms{phantomIdx}.name, [1:length(obj.phantoms)], 'UniformOutput', false);

            if isa(phantom, 'baseData_genericPhantom')
                if ~any(strcmp(currentPhantomNames, phantom.name))
                    obj.phantoms = [obj.phantoms, {phantom}];
                    obj.nPhantoms = obj.nPhantoms +1;
                else
                    matRad_cfg.dispWarning('Phantom already present');
                end
            else
                matRad_cfg.dispError('Given class is not a phantom');
            end
        end

        function removePhantom(obj,phantom)
            matRad_cfg = MatRad_Config.instance();
            currentPhantomNames = arrayfun(@(phantomIdx) obj.phantoms{phantomIdx}.name, [1:length(obj.phantoms)], 'UniformOutput', false);
            if ~isempty(obj.phantoms)
                 if any(strcmp(phantom.name, currentPhantomNames))
                    phantomIdx = find(strcmp(phantom.name, currentPhantomNames));
                    obj.phantoms(phantomIdx) = [];
                    obj.nPhantoms = obj.nPhantoms -1;
                else
                    matRad_cfg.dispError(['No phantom named: ', phantomName, ' found']);
                end
            else
                matRad_cfg.dispError('No phantom to be removed');
            end
        end
    end

end