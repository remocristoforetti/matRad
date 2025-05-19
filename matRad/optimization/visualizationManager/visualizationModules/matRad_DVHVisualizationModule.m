classdef matRad_DVHVisualizationModule < matRad_VisualizationModule

    properties
        quantity;
        structures;
        cst;
        
        dataRequest = {'d'};
        data;
    end

    methods
        function this = matRad_DVHVisualizationModule(quantity, cst, structures)
            if ~exist('structures', 'var')
                structures = [1:size(cst,1)];
            end


            this.quantity = quantity;
            this.cst = cst;
            this.structures = structures;
        end

        function updateData(this, dataInput)
       
            this.data.d = dataInput.d.(this.quantity){1};            
        
        end

        function isValid = validate(this)

            % Validate the model
            isValid = true;

        end

        function plotData(this)

            dvh = matRad_calcDVH(this.cst,this.data.d);

            hold(this.hAx, 'off');
            for i=this.structures
                plot(this.hAx, dvh(i).doseGrid, dvh(i).volumePoints,'.-', 'Color', this.cst{i,5}.visibleColor, 'DisplayName',this.cst{i,3});
                hold(this.hAx, 'on');
            end
            grid(this.hAx, 'on');
            legend(this.hAx);
            xlabel(this.hAx,'Dose [Gy]');
            ylabel(this.hAx, 'Volume [%]');
    
        end
    end
end
