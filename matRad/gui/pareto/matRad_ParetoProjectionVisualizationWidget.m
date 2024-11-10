classdef matRad_ParetoProjectionVisualizationWidget < matRad_Widget
    
    properties
        %hFig;
    end

    methods
        function this = matRad_ParetoProjectionVisualizationWidget(handleParent)
            this = this@matRad_Widget(handleParent);
        end

        function h = initialize(this)

            %hFig, fInd, cst, optiProb, additionalPoint
            % if isempty(this.hFig)
            %     this.hFig = %axes(this.widgetHandle,'Position',[0 0 1 1]);
            % end
            currentPoint = evalin('base','ParetoHelperObject.currentPoint');
            cst = evalin('base','cst');
            optiProb = evalin('base', 'retStruct.optiProb');
            fInds = evalin('base', 'retStruct.finds');
            matRad_projectParetoSurface(this.widgetHandle,fInds, cst, optiProb, currentPoint);

        end

    end

    methods (Access = protected)
        function this = createLayout(this)
            this.createHandles();
        end
        
    end
end