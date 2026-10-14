classdef matRad_StoppingPowerModel < matRad_TabulatedDoseAveragedModel
% This class dose-averages tabulated, fragment- and energy-resolved
% alpha and beta values (e.g. from LEM) using a stopping power table for
% weighting.
%
% quantityTable holds one entry per available cell line: either an
% already loaded table (struct with fields meta, data) or the name of a
% table to be lazily loaded from the alphaBetaTables folder on first use
% (see loadAlphaBetaTable). Each table's meta.modelParameters.cellLine
% identifies the cell line it was fitted for, and is used by
% getQuantityTable to match a table against a requested cell line.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2023-2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    properties (Constant)
        model                  = 'doseAveragedStoppingPower'
        requiredQuantities     = {}
        matrixContainer        = {'mSPDose'};
    end

    properties
        quantityName           = {'dEdx'}
        quantityTable
    
        storeQuantityIJ        = false
    end

    methods


        function this = matRad_StoppingPowerModel()
            this@matRad_TabulatedDoseAveragedModel()

            this.quantityTable = {};
        end

        function bixel = calcBiologicalQuantitiesForBixel(this,bixel,kernel)
            bixel = calcBiologicalQuantitiesForBixel@matRad_LQBasedModel(this,bixel);
            
            for  cellLineIdx = 1:numel(unique(bixel.cellLineIndex))
                voxelIdx = find(bixel.cellLineIndex == cellLineIdx);
                bixel.alpha(voxelIdx) = kernel.alpha(voxelIdx,cellLineIdx);
                bixel.beta(voxelIdx)  = kernel.beta(voxelIdx,cellLineIdx);
            
            end

            if this.storeQuantityIJ
                bixel.(this.matrixContainer{1}) = kernel.dEdx .* bixel.physicalDose;
            end

        end


        function kernel = computeDoseAveragedKernel(this, machine, cellLines)

            kernel = computeDoseAveragedKernel@matRad_TabulatedDoseAveragedModel(this, machine, cellLines);

            kernel.alpha = NaN*ones(size(kernel.dEdx,1), numel(cellLines));
            kernel.beta  = NaN*ones(size(kernel.dEdx,1), numel(cellLines));

        end

        function tables = getQuantityTable(this,~)
            % Resolves and returns the alpha/beta tables matching the
            % given cell line(s), using table.meta.modelParameters.cellLine,
            % in the same order as cellLines. Warns for any requested cell
            % line that isn't found.
            %
            % As a side effect, this.cellLines is set to the subset of the
            % requested cellLines that were actually found in the table,
            % preserving the order they were passed in.
            matRad_cfg = MatRad_Config.instance();

            tables = {this.stoppingPowerTable};
 
        end

        function validateQuantityTable(this)

        end

    end

    methods (Static)

    end
end