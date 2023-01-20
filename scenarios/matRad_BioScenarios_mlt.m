classdef matRad_BioScenarios_mlt < matRad_BioScenarios
%  matRad_RandomScenarios
%  Implements randomly sampled scenarios
%
% constructor
%   matRad_RandomScenarios()
%   matRad_RandomScenarios(ct)
%
% input
%   ct:                 ct cube
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2022 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
      SelectedBioModels;
      nSamples = 1;
    end

    properties (SetAccess=protected)
       name = 'BioScenarios_MultiModels';
    end
    methods
        function this = matRad_BioScenarios_mlt(ct)
           if nargin == 0
                superclassArgs = {};
            else
                superclassArgs = {ct};
            end
            this@matRad_BioScenarios(superclassArgs{:});
            this.updateScenarios();
        end

        %% Setters & Update
        
        function set.SelectedBioModels(this,value)
           matRad_cfg = MatRad_Config.instance();
           nModels = size(value,2);
           SelBioModels = [];


           for k=1:nModels
              if any(strcmp(value(k),this.AvailableBioModels))
                  SelBioModels = [SelBioModels,value(k)];
              else
                 matRad_cfg.dispWarning('Model %s not available', value{k});
              end
           end
           this.SelectedBioModels = SelBioModels;
           if size(SelBioModels,2) >1
               this.nSamples = size(SelBioModels,2);
           else %BioModel might not be defined
              this.ParameterDeviations = zeros(this.BioModel.NumModelParameters);
              %why is this here? Should not care of ParameterDeviations if
              %multiModel
           end
        end
        
        function set.nSamples(this,nSamples)
            valid = isnumeric(nSamples) && isscalar(nSamples) && mod(nSamples,1) == 0 && nSamples > 0;
            matRad_cfg = MatRad_Config.instance();
            if ~valid 
                matRad_cfg.dispError('Invalid value for nSamples! Needs to be a positive integer!');
            end
            
            if ~isempty(this.SelectedBioModels)
               if size(this.SelectedBioModels,2) == nSamples
                  this.nSamples = nSamples;
               else
                  matRad_cfg.dispWarning('nSamples should match the number of given models');
                  this.nSamples = size(this.SelectedBioModels,2);
               end
            else
               this.nSamples = nSamples;
            end
            this.updateScenarios();
        end
        
        function Tp = calcTissueParameter(this,cst,numOfVoxels,ctScen)
               for k =1:this.nSamples
                   BioModel = this.BioModels{k};
                   Tp(k) = {BioModel.calcTissueParameter(cst,numOfVoxels,ctScen)};
               end    

        end
        
        
        function scenarios = updateScenarios(this)
            matRad_cfg = MatRad_Config.instance();
            
            %Multivariate Normal Sampling
            Sigma = diag([this.shiftSD,this.rangeAbsSD,this.rangeRelSD./100].^2);
            d = size(Sigma,1);
            [cs,p] = chol(Sigma);
            

            scenarios = zeros(this.nSamples,5);
            
            %Scenario Probability from pdf
            this.scenProb = (2*pi)^(-d/2) * exp(-0.5*sum((scenarios/cs).^2, 2)) / prod(diag(cs));
            this.scenWeight = ones(this.nSamples,1)./this.nSamples;
            
            this.scenForProb = scenarios;
            %set variables
            this.totNumShiftScen = this.nSamples;
            this.totNumRangeScen = this.nSamples;
            this.totNumScen = this.nSamples; %check because of CT scenarios
            
            %Individual shifts
            this.relRangeShift = scenarios(:,5);
            this.absRangeShift = scenarios(:,4);
            this.isoShift = scenarios(:,1:3);

            this.maxAbsRangeShift = max(this.absRangeShift);
            this.maxRelRangeShift = max(this.relRangeShift);

            %Mask for scenario selection
            this.scenMask = false(this.numOfCtScen,this.totNumShiftScen,this.totNumRangeScen);

            for sCt = 1:this.numOfCtScen
                this.scenMask(sCt,:,:) = diag(true(this.nSamples,1));
            end
            
            [x{1}, x{2}, x{3}] = ind2sub(size(this.scenMask),find(this.scenMask));
            this.linearMask    = cell2mat(x);
            totNumScen    = sum(this.scenMask(:));

            if totNumScen ~= this.totNumScen
                matRad_cfg.dispWarning('Check Implementation of Total Scenario computation - given %d but found %d!',this.totNumScen,totNumScen);
                this.totNumScen = totNumScen;
            end
            
            this.updateBioModels(this.SelectedBioModels);
        end
        
        
        function this = updateBioModels(this,BioModels)
           nBioModels = size(BioModels,2);
           BioModelInstances = [];
           for k=1:nBioModels
               BioModelInstances = [BioModelInstances, this.setBioModel(BioModels{k})];
           end
           this.BioModels = BioModelInstances;
        end
    end


end