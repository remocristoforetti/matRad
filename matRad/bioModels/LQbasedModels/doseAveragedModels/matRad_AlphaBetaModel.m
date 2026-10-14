classdef matRad_AlphaBetaModel < matRad_TabulatedDoseAveragedModel
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
        model                  = 'doseAveragedAlphaBetaModel'
        requiredQuantities     = {}                                      % kernels in base data needed for the alpha/beta calculation
    end

    properties
        cellLines
        quantityTable % cell array, one entry per cell line: either a loaded table (struct with fields meta, data) or the name of a table to load
    end

    properties (SetAccess= protected)
        storeQuantityIJ        = false;
    end

    properties
        quantityName           = {'alpha'; 'sqrtBeta'}
    end

    methods


        function this = matRad_AlphaBetaModel()
            this@matRad_TabulatedDoseAveragedModel()

            this.cellLines = [];
            this.quantityTable = {};
        end

        function bixel = calcBiologicalQuantitiesForBixel(this,bixel,kernel)
            bixel = calcBiologicalQuantitiesForBixel@matRad_LQBasedModel(this,bixel);
            
            for cellLineIdx = 1:numel(this.cellLines)
                voxelIdx = find(bixel.cellLineIndex == cellLineIdx);
                bixel.alpha(voxelIdx) = kernel.alpha(voxelIdx,cellLineIdx);
                bixel.beta(voxelIdx)  = kernel.sqrtBeta(voxelIdx,cellLineIdx).^2;
            end

            % Voxels outside dose grid 
            bixel.alpha(isnan(bixel.alpha)) = 0;
            bixel.beta(isnan(bixel.beta)) = 0;
        end

        function kernel = computeDoseAveragedKernel(this, machine, cellLines)
            kernel = computeDoseAveragedKernel@matRad_TabulatedDoseAveragedModel(this, machine, cellLines);

            kernel.beta = kernel.sqrtBeta.^2;
        end

        function tables = getQuantityTable(this, cellLines)
            % Resolves and returns the alpha/beta tables matching the
            % given cell line(s), using table.meta.modelParameters.cellLine,
            % in the same order as cellLines. Warns for any requested cell
            % line that isn't found.
            %
            % As a side effect, this.cellLines is set to the subset of the
            % requested cellLines that were actually found in the table,
            % preserving the order they were passed in.
            matRad_cfg = MatRad_Config.instance();

            this.validateQuantityTable();

            availableCellLines = cellfun(@(t) char(t.meta.modelParameters.cellLine), this.quantityTable, 'UniformOutput', false);

            cellLineNames = {cellLines(:).cellLineName};
            isAvailable = ismember(cellLineNames,availableCellLines);

            this.cellLines = cellLines(isAvailable);

            tables = cell(size(cellLines));
            for c = 1:numel(cellLines)
                if isAvailable(c)
                    matchIx   = find(strcmp(availableCellLines,cellLineNames),1,'first');
                    tables{c} = this.quantityTable{matchIx};
                    for j=1:numel(tables{c}.data)
                        tables{c}.data(j).sqrtBeta = sqrt(tables{c}.data(j).beta);
                    end
                else
                    matRad_cfg.dispWarning('Cell line ''%s'' was requested but not found!',cellLines{c});
                end
            end
        end

        function validateQuantityTable(this)
            % Checks that quantityTable is set, loads any not-yet-loaded
            % entry (given as a filename) via loadAlphaBetaTable, and
            % sanity-checks the resulting table structure
            matRad_cfg = MatRad_Config.instance();

            validateQuantityTable@matRad_TabulatedDoseAveragedModel(this);

            for t = 1:numel(this.quantityTable)

                % This part needs to be included, defines how to load the
                % table
                if ischar(this.quantityTable{t}) || isstring(this.quantityTable{t})
                    this.quantityTable{t} = matRad_AlphaBetaModel.loadAlphaBetaTable(char(this.quantityTable{t}));
                elseif ~isstruct(this.quantityTable{t})
                    matRad_cfg.dispError('quantityTable entry %d must be either a loaded table struct or the name of a table to load!',t);
                end

                if ~isfield(this.quantityTable{t},'meta') || ~isfield(this.quantityTable{t},'data')
                    matRad_cfg.dispError('quantityTable entry %d does not contain the required meta/data fields!',t);
                end

                if ~isfield(this.quantityTable{t}.meta,'modelParameters') || ~isfield(this.quantityTable{t}.meta.modelParameters,'cellLine')
                    matRad_cfg.dispError('quantityTable entry %d is missing meta.modelParameters.cellLine!',t);
                end

                tableFragments = [ [this.quantityTable{t}.data.Z]', [this.quantityTable{t}.data.A]' ];
                isAvailable = ismember([this.selectedFragments.Z(:), this.selectedFragments.A(:)], tableFragments, 'rows');

                if ~all(isAvailable)
                    missingFragments = arrayfun(@(z,a) sprintf('Z=%d,A=%d',z,a), this.selectedFragments.Z(~isAvailable), this.selectedFragments.A(~isAvailable), 'UniformOutput', false);
                    matRad_cfg.dispError('The following selected fragments are not available in quantityTable entry %d: %s',t,strjoin(missingFragments,'; '));
                end
            end
        end

    end

    methods (Static)

        function table = loadAlphaBetaTable(fileName)
            % Loads and returns the table with name: fileName
            matRad_cfg = MatRad_Config.instance();

            searchPath = {fullfile(matRad_cfg.matRadSrcRoot,'bioModels','alphaBetaTables'), ... % default matRad folder
                          fullfile(matRad_cfg.primaryUserFolder,'alphaBetaTables')};             % user defined table

            [inputPath, inputFileName, inputExtension] = fileparts(fileName);
            
            table = [];
            if isempty(inputPath)
                for p = 1:numel(searchPath)
                    try
                        tmp = load(fullfile(searchPath{p}, [fileName, '.mat']),'alphaBetaTable');
                        table = tmp.alphaBetaTable;
                        break;
                    catch
                        continue;
                    end
                end
            else
                try
                    tmp = load(fullfile(inputPath, sprintf('%s%s', inputFileName, inputExtension)));
                    table = tmp.alphaBetaTable;
                catch
                end
            end

            if isempty(table)
                matRad_cfg.dispError('Cannot find alpha/beta table: %s',fileName);
            end

            if ~isstruct(table) || ~isequal(sort(fieldnames(table)),sort({'meta';'data'}))
                matRad_cfg.dispError('Provided alpha/beta table does not contain the required meta/data fields!');
            end

            matRad_cfg.dispInfo('Loaded alpha/beta table ''%s'' with %d fragment(s).\n',fileName,numel(table.data));
        end

    end
end