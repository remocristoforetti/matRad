classdef matRad_DistributionVisModule < matRad_VisualizationModule

    properties
        distributionParameters;

        
        dataRequest = {'d'};
        data;
    end

    methods
        function this = matRad_DistributionVisModule(distributionParameters)
            if exist('distributionParameters', 'var')
                this.distributionParameters = distributionParameters;
            end
        end

        function updateData(this, dataInput)
        
            
            this.data.d = dataInput.d;
        
        end

        function isValid = validate(this)

            % Validate the model
            isValid = true;

        end

        function plotData(this)

            doseGrid = this.distributionParameters.doseGrid;
            ctGrid   = this.distributionParameters.ct;
            
            dDoseGrid = reshape(this.data.d.(this.distributionParameters.quantity){1}, doseGrid.dimensions);
            dCtGrid = matRad_interp3(doseGrid.x,doseGrid.y',doseGrid.z, ...
                                     dDoseGrid, ...
                                     ctGrid.x,ctGrid.y',ctGrid.z,'linear',0);
             
            matRad_plotSliceWrapper(this.hAx, this.distributionParameters.ct, this.distributionParameters.cst, 1, dCtGrid, this.distributionParameters.plane, this.distributionParameters.slice, [],[],[],[],[],[],ones(1,size(this.distributionParameters.cst,1)));
    
        end
    end
end
