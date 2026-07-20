classdef (Abstract) matRad_TabulatedDoseAveragedModel < matRad_LQBasedModel
% This is an Abstract class implementing dose-averaging of a tabulated,
% fragment- and energy-resolved quantity (e.g. LET, alpha, ...) using a
% tabulated stopping power table for weighting.
%
% Concrete subclasses must implement:
%   getQuantityTable(this,loopVar)  returns the table matching the given
%                               loop variable (e.g. a cell line name, or
%                               whatever grouping is meaningful to the
%                               subclass; opaque to this class). The
%                               table is a struct array, one entry per
%                               fragment, with fields Z, A, energies, and
%                               a field named after quantityName holding
%                               the values (vs. energies) to be
%                               dose-averaged. May return empty if no
%                               table matches the given loop variable.
% and declare the constant property:
%   quantityName                name of the field (in the struct array
%                               returned by getQuantityTable) holding the
%                               quantity to average. Also used as the
%                               field name of the resulting kernel.
%
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

    properties (Abstract, Constant)
        quantityName           % name of the averaged quantity / field (see class header)
    end

    properties (Abstract)
        quantityTable
    end

    properties (Constant)
        possibleRadiationModes = {'protons', 'carbon', 'helium', 'oxigen'}  % radiation modalitites compatible with the model
    end

    properties
        stoppingPowerTable                 % struct with fields: meta, data (array of structs with Z,A,energy,dEdx)
    
        selectedFragments    % indices into SPtable.data of fragments included for averaging
    end

    properties (SetAccess = protected)
        interpolatedStoppingPower
    end

    methods (Abstract)
        % Returns the table matching the given loop variable (e.g. a cell
        % line name, or whatever grouping is meaningful to the subclass;
        % opaque to this class), or empty if none matches. The table is
        % a struct array, one entry per fragment, with fields Z, A,
        % energies, and this.quantityName (or one field per name, if
        % quantityName lists several quantities), defining the
        % energy-resolved quantity/quantities to be dose-averaged.
        % Implemented by concrete subclasses, typically by matching the
        % loop variable against this.quantityTable (loading any
        % not-yet-loaded entries first).
        quantityTable = getQuantityTable(this, loopVariable)
    end

    methods
        function this = matRad_TabulatedDoseAveragedModel
            this@matRad_LQBasedModel()
            this.stoppingPowerTable = [];
            this.selectedFragments = [];
        end

        function loadSPTable(this, filename)
            % Load Stopping power table in fileName
            matRad_cfg = MatRad_Config.instance();

            searchPath = {fullfile(matRad_cfg.matRadSrcRoot,'bioModels','StoppingPowerTables'), ... % default matRad folder
                          fullfile(matRad_cfg.primaryUserFolder,'StoppingPowerTables')};             % user defined table

            loadedTable = [];
            for p = 1:numel(searchPath)
                try
                    tmp = load(fullfile(searchPath{p}, [filename, '.mat']),'SPtable');
                    loadedTable = tmp.SPtable;
                    break;
                catch
                    continue;
                end
            end

            if isempty(loadedTable)
                matRad_cfg.dispError('Cannot find stopping power table: %s',filename);
            end

            if ~isstruct(loadedTable) || ~isequal(sort(fieldnames(loadedTable)),sort({'meta';'data'}))
                matRad_cfg.dispError('Provided stopping power table does not contain the required meta/data fields!');
            end

            this.stoppingPowerTable = loadedTable;

            matRad_cfg.dispInfo('Loaded stopping power table ''%s'' with %d fragment(s).\n',filename,numel(loadedTable.data));
        end

        function data = selectFragmentsFromTable(this,table, Zlist, Alist)
            % Returns the indices into table.data of the entries matching
            % the given (Z,A) fragment pairs (paired element-wise).
            % Fragments that are not found are returned as NaN.
            matRad_cfg = MatRad_Config.instance();

            if numel(Zlist) ~= numel(Alist)
                matRad_cfg.dispError('Zlist and Alist must have the same number of elements!');
            end

            tableZ = [table.data.Z];
            tableA = [table.data.A];

            data = NaN(numel(Zlist),1);
            for k = 1:numel(Zlist)
                ix = find(tableZ == Zlist(k) & tableA == Alist(k),1,'first');
                if ~isempty(ix)
                    data(k) = ix;
                end
            end
        end

        function interpolateFragmentSP(this, energies, Z,A)
            % Interpolate SPtable for the fragment with Z,A, at the specified energies.
            % - select fragment from this.SPtable
            % - interpolate the dEdx
            % - save in interpolatedStoppingPower for the specific fragment. interpolatedStoppingPower has the same structure as SPtable but only for the selectedFragments
            matRad_cfg = MatRad_Config.instance();

            ix = this.selectFragmentsFromTable(this.stoppingPowerTable,Z,A);

            if isnan(ix)
                matRad_cfg.dispError('Fragment Z=%d, A=%d is not available in the loaded stopping power table!',Z,A);
            end

            fragmentData = this.stoppingPowerTable.data(ix);

            newEntry.Z        = Z;
            newEntry.A        = A;
            newEntry.energies = energies(:);
            newEntry.dEdx     = matRad_interp1(fragmentData.energies(:), fragmentData.dEdx(:), energies(:), 'extrap');

            if isempty(this.interpolatedStoppingPower)
                this.interpolatedStoppingPower = newEntry;
            else
                existingIx = find([this.interpolatedStoppingPower.Z] == Z & [this.interpolatedStoppingPower.A] == A,1,'first');
                if ~isempty(existingIx)
                    this.interpolatedStoppingPower(existingIx) = newEntry;
                else
                    this.interpolatedStoppingPower(end+1) = newEntry;
                end
            end
        end

        function validateStoppingPowerTable(this)
            % Check if the selectedFragments are available in the SPtable
            matRad_cfg = MatRad_Config.instance();

            if isempty(this.stoppingPowerTable)
                matRad_cfg.dispError('No stopping power table loaded!');
            elseif ischar(this.stoppingPowerTable)
                this.loadSPTable(this.stoppingPowerTable);
            end

            if isempty(this.selectedFragments) || ~isfield(this.selectedFragments,'Z') || ~isfield(this.selectedFragments,'A') ...
                    || isempty(this.selectedFragments.Z) || isempty(this.selectedFragments.A)
                matRad_cfg.dispError('No fragments selected for dose-averaging! selectedFragments must be a struct with fields Z and A.');
            end

            if numel(this.selectedFragments.Z) ~= numel(this.selectedFragments.A)
                matRad_cfg.dispError('selectedFragments.Z and selectedFragments.A must have the same number of elements!');
            end

            tableFragments = [ [this.stoppingPowerTable.data.Z]', [this.stoppingPowerTable.data.A]' ];
            isAvailable = ismember([this.selectedFragments.Z(:), this.selectedFragments.A(:)], tableFragments, 'rows');

            if ~all(isAvailable)
                missingFragments = arrayfun(@(z,a) sprintf('Z=%d,A=%d',z,a), this.selectedFragments.Z(~isAvailable), this.selectedFragments.A(~isAvailable), 'UniformOutput', false);
                matRad_cfg.dispError('The following selected fragments are not available in the loaded stopping power table: %s',strjoin(missingFragments,'; '));
            end

        end

        function validateQuantityTable(this)
            % Checks that quantityTable is set, loads any not-yet-loaded
            % entry (given as a filename) via loadAlphaBetaTable, and
            % sanity-checks the resulting table structure
            matRad_cfg = MatRad_Config.instance();
    
            if ~iscell(this.quantityTable)
                this.quantityTable = {this.quantityTable};
            end
    
            if isempty(this.quantityTable)
                matRad_cfg.dispError('No quantity table(s) set for dose-averaging!');
            end
    
            if isempty(this.selectedFragments) || ~isfield(this.selectedFragments,'Z') || ~isfield(this.selectedFragments,'A') ...
                    || isempty(this.selectedFragments.Z) || isempty(this.selectedFragments.A)
                matRad_cfg.dispError('No fragments selected for dose-averaging! selectedFragments must be a struct with fields Z and A.');
            end
    
            if numel(this.selectedFragments.Z) ~= numel(this.selectedFragments.A)
                matRad_cfg.dispError('selectedFragments.Z and selectedFragments.A must have the same number of elements!');
            end

        end

        function kernel = computeDoseAveragedKernel(this, machine, cellLines)
            % The machine input is a structure with fields:
            % spectra: a struct array, one element for each fragment with fields:
                % - Z
                % - A
                % - energies
                % - fluence: 2D matrix of dimension: number of energies vs number of depths
            %
            % cellLineNames is only used here as a loop variable: for each entry, the
            % matching quantity table is requested from the subclass via getQuantityTable(cellLine).
            % This class attaches no further meaning to it.

            % This function does:

            % Sanity check on the input machine and the stopping power table
            % Gets, for each requested loop variable (cellLineNames), the quantity table to
            % average (an energy vs quantity table) defined by specific subclasses of this class
            % interpolates the stopping power table (SPtable) on the energies provided in the machine for the required fragments
            % computes the product stopping power times fluence (in the machine) times the quantity (summing over all the energies and fragments), get an array that has the dimension of number of depths
            % computes the product stopping power times fluence (in the machine) (summing over all the energies and fragments),get an array that has the dimension of number of depths
            % takes the ratio of the previous two arrays, separately for each requested loop variable
            matRad_cfg = MatRad_Config.instance();

            % Sanity checks
            if ~isfield(machine,'fluenceSpectra') || isempty(machine.fluenceSpectra)
                matRad_cfg.dispError('Machine base data does not contain fragment spectra required for dose-averaging!');
            end

            this.validateStoppingPowerTable();

            quantityNames = cellstr(this.quantityName); % quantityName may be a single name or a cell array of names

            % resolve the quantity table for each requested loop variable once upfront; a
            % subclass may return empty if no table matches a given entry
            % quantityTables = cell(1,numel(cellLineNames));
            % for c = 1:numel(cellLineNames)
            quantityTables = this.getQuantityTable(cellLines);
            % end

            machineZ = [machine.fluenceSpectra.Z];
            machineA = [machine.fluenceSpectra.A];

            numOfDepths = size(machine.fluenceSpectra(1).fluence,1);
            numOfTables = numel(quantityTables); % num of cellLines

            numerator   = zeros(numOfDepths,numel(quantityNames),numOfTables);
            denominator = zeros(numOfDepths,1);

            for idx = 1:numel(this.selectedFragments.Z)
                Z = this.selectedFragments.Z(idx);
                A = this.selectedFragments.A(idx);

                fragmentIx = find(machineZ == Z & machineA == A,1,'first');
                if isempty(fragmentIx)
                    matRad_cfg.dispError('Selected fragment Z=%d, A=%d for dose-averaging not found in the provided base data spectra!',Z,A);
                end
                fragmentSpectrum = machine.fluenceSpectra(fragmentIx);

                % interpolate the stopping power on the fragment's fluence energy grid (shared across tables)
                this.interpolateFragmentSP(fragmentSpectrum.energyBins,Z,A);
                spIx = find([this.interpolatedStoppingPower.Z] == Z & [this.interpolatedStoppingPower.A] == A,1,'first');
                sp = this.interpolatedStoppingPower(spIx).dEdx(:);

                weight = sp' .*fragmentSpectrum.fluence; % nEnergies x nDepths
                denominator = denominator + sum(weight,2);

                for t = 1:numOfTables
                    currentTable = quantityTables{t}.data;
                    if isempty(currentTable)
                        continue; % no table found for this loop variable, subclass already warned
                    end

                    if ~all(isfield(currentTable,[{'Z','A','energies'}, quantityNames]))
                        matRad_cfg.dispError('Quantity table for ''%s'' provided by %s does not contain the required fields Z, A, energies, %s!',cellLines{t},class(this),strjoin(quantityNames,', '));
                    end

                    quantityIx = find([currentTable.Z] == Z & [currentTable.A] == A,1,'first');
                    if isempty(quantityIx)
                        matRad_cfg.dispError('Selected fragment Z=%d, A=%d not found in the quantity table for ''%s''!',Z,A,cellLines{t});
                    end
                    fragmentQuantity = currentTable(quantityIx);

                    for q = 1:numel(quantityNames)
                        % interpolate the quantity to average on the same energy grid, in case it was tabulated on a different grid
                        quantity = matRad_interp1(fragmentQuantity.energies(:), fragmentQuantity.(quantityNames{q})(:), fragmentSpectrum.energyBins(:), 'extrap');
                        numerator(:,q,t) = numerator(:,q,t) + weight * quantity;
                    end
                end
            end

            for q = 1:numel(quantityNames)
                % depths x loop variables, in the same order as cellLineNames
                kernel.(quantityNames{q}) = squeeze(numerator(:,q,:)) ./ denominator(:);
            end
        end

    end

    methods (Static)
        % Static helpers can be added here in future
    end
end