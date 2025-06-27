classdef matRad_TabulatedSpectralKernelBasedModel < matRad_LQRBETabulatedModel
% This is class implementig a spectra-based tabulated RBE model.
% Spectral kernels should be provided in the base data
% The model can handle multiple tissue alphaX/betaX ratio specified by the 
% cst structure, as long as a compatible RBEtable is provided.
%
% Properties of the model that can be defined by the users:
%   RBEtableName        filename of the table to be included
%
%   fragmentsToInclude  specifies which fragments to include in the
%                       biological calculation (i.e. 'H', 'He', 'C',...)
%
%   weightBy            specifies which kernel data is used in combination
%                       with the RBE table (i.e. 'Fluence', 'EnergyDeposit')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2023 the matRad development team.
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
        model = 'TAB';
        
        requiredQuantities = {'spectra','physicalDose'};
        kernelQuantities = {'spectra'};
        possibleRadiationModes = {'protons','helium','carbon'};
    end

    properties
        weightBy;
        baseDataKernel; % This is set upon initialization of the model according to the specific weighting quantity
    end

    methods
        function this = matRad_TabulatedSpectralKernelBasedModel()
            this@matRad_LQRBETabulatedModel();
            this.assignDefaultProperties();
        end


        function bixel = calcBiologicalQuantitiesForBixel(this,bixel,kernel)
        
            % This function implements the variable alpha/beta calculation
            % for the bixel.
            % These values are linearly combined for each fragment and
            % energy. The weighting factor is specified by the selected
            % spectra kernel.
            % Formulation for the mixed field:
            % alpha{d} = sum(sum(alpha(E,F).*Phi(E,F,d))./sum(sum(Phi(E,F,d)));
            % with:
            %   d           radiological depth
            %   E           energy
            %   F           fragment (particle class)
            %   Phi(E,F,d)  Spectum (i.e. Fluence)
            %   alpha(E,F)  alpha value from table for energy E and frgmant F
            % sums are performed over energy and fragments


            bixel = calcBiologicalQuantitiesForBixel@matRad_LQRBETabulatedModel(this,bixel);

            nFragments = numel(this.fragmentsToInclude);

            % Collect the spectra from the bixel, for each fragment
            % For now assume all fragments have same energy.
            spectraEnergies = bixel.baseData.spectra.(this.weightBy).energies;
            
            % Get the tissue classes within the bixel
            bixelTissueIndexes = unique(bixel.vTissueIndex)';
            
            % Interpolate the alpha/beta table for the specific fragment (including primaries)
            % and alphaX/betaX ratio (tissue class)
            % Change here how this is saved into alphaE to have separate
            % fragment fields
            for i=bixelTissueIndexes
                [curTissueAlphaE(i,:), curTissueBetaE(i,:)] = cellfun(@(fragment) this.interpolateRBETableForBixel(spectraEnergies, fragment, i),this.fragmentsToInclude, 'UniformOutput',false);          
            end

            alphaE = arrayfun(@(fragment) horzcat(curTissueAlphaE{:,fragment}), [1:numel(this.fragmentsToInclude)], 'UniformOutput',false);
            betaE  = arrayfun(@(fragment) horzcat(curTissueBetaE{:,fragment}),  [1:numel(this.fragmentsToInclude)], 'UniformOutput',false);

            alphaE = cell2struct(alphaE, this.fragmentsToInclude, 2);
            betaE  = cell2struct(betaE,  this.fragmentsToInclude, 2);

            % Get the spectra kernels to be used (one for each fragment).
            bixelSpectra = cell2struct(cellfun(@(fragment) squeeze(kernel.(fragment)),this.fragmentsToInclude, 'UniformOutput',false), this.fragmentsToInclude, 2);
           
            % Get normalization for each fragment. This is sum over
            % energies
            spectraFragmentDenominator = cellfun(@(fragment) sum(bixelSpectra.(fragment),2),this.fragmentsToInclude, 'UniformOutput',false);

            % Get total normalization by summing over the fragments
            spectraBixelDenominator = sum([spectraFragmentDenominator{:}],2);

            % Create containers for alpha and beta for each fragment
            alphaWeightedFragmentSpectra = NaN*ones(numel(bixel.radDepths),nFragments);
            betaWeightedFragmentSpectra  = NaN*ones(numel(bixel.radDepths),nFragments);
            
            
            for i=bixelTissueIndexes
 
                % Get the spectrum-weighted alpha/beta for each fragment (sums over the energy bins at each bixel radDepth)
                tmpAlphaWeightedFragmentSpectra = cellfun(@(fragment) bixelSpectra.(fragment)*alphaE.(fragment)(:,i), this.fragmentsToInclude, 'UniformOutput',false);
                tmpBetaWeightedFragmentSpectra  = cellfun(@(fragment) bixelSpectra.(fragment)*betaE.(fragment)(:,i), this.fragmentsToInclude, 'UniformOutput',false);
                
                % Put it in matrix form and sum over the fragments. Only
                % include values for the correct tissue class for each
                % radDepth
                matTmpAlphaWeightedFragmentSpectra = [tmpAlphaWeightedFragmentSpectra{:}];
                matTmpBetaWeightedFragmentSpectra  = [tmpBetaWeightedFragmentSpectra{:}];

                alphaWeightedSpectra(bixel.vTissueIndex == i,:) = sum(matTmpAlphaWeightedFragmentSpectra(bixel.vTissueIndex == i,:),2);
                betaWeightedSpectra(bixel.vTissueIndex == i,:)  = sum(matTmpBetaWeightedFragmentSpectra(bixel.vTissueIndex == i,:),2);
            end

            % Normalize the quantities
            bixel.alpha = alphaWeightedSpectra./spectraBixelDenominator;
            bixel.beta  = betaWeightedSpectra./spectraBixelDenominator;
        end
    end

    methods

        function assignDefaultProperties(this)
            
            assignDefaultProperties@matRad_LQRBETabulatedModel(this);
            this.weightBy = 'Fluence';
            
        end


        function set.weightBy(this, value)

            if ischar(value)
                this.weightBy = value;
                this.updatePropertyValues();
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Unrecognized weighting quantity: %s', tostring(value));
            end

        end

       
        function updatePropertyValues(this)
        
            if ~isempty(this.weightBy) && ~isempty(this.fragmentsToInclude)
                
                this.baseDataKernel = cellfun(@(fragment) ['spectra.', this.weightBy, '.', fragment, '.Data'], this.fragmentsToInclude, 'UniformOutput',false); %[requiredSpectraData;requriedEnergies];
            end


            if ~isempty(this.RBEtable)
                this.checkTableConsistency(this.RBEtable, this.fragmentsToInclude);
            end

        end

    end
end