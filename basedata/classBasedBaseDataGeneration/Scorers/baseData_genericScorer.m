classdef baseData_genericScorer < handle
    properties
        name;
        binnedInEnergy;
        fragments = {'none'};
        energyBinning; %a class for every fragment

    end

    properties (SetAccess=protected)

    end

    methods
        function obj = baseData_genericScorer()

        end


        function addFragment(obj, fragment)
            if isa(fragment, 'string')
                switch fragment

                    case 'proton'
                    obj.fragments = [obj.fragments, baseData_ion_proton()];
                    case 'none'

                end
            elseif isa(fragment, 'baseData_ion')
                obj.fragments = [obj.fragments, fragment];
            else

            end
        end
    end
end