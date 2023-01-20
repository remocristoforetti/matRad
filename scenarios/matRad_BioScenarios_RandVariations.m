classdef matRad_BioScenarios_RandVariations < matRad_BioScenarios
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
      ParameterStandardDeviations;
      ModelParameters;
      %name;% = 'bioScen_rndV';
      includeNominalScenario = true;
    end
    
    properties (SetAccess=protected)
       name = 'BioScenarios_RandVariations';
    end

    properties (Hidden)
      NumModelParameters;
    end
    methods
        function this = matRad_BioScenarios_RandVariations(ct)
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
           if ~strcmp(value,this.AvailableBioModels)
               matRad_cfg.dispWarning('Model %s not available', value{k});
           else
              this.SelectedBioModels = value;
           end
           
           this.updateScenarios();
        end

        function Variations = calcParameterVariations(this,TissueParameters)
               Fields = this.ModelParameters;%fieldnames(TissueParameters);
               
               nFields = this.NumModelParameters; %sum(~strcmp(Fields,'tissueClass'));
               %! Fields and fields of tissueParameters might not match
               nTissues = size(getfield(TissueParameters,Fields{1}),2);
               Variations = TissueParameters;
               for k=1:nFields

                  if ~strcmp(Fields{k},'tissueClass')
                     RandomVariation = normrnd(getfield(TissueParameters,Fields{k}),this.ParameterStandardDeviations(k)*ones(1,nTissues));
                     RandomVariation(RandomVariation<0)=0;
                     %Variations = setfield(Variations,Fields{k},[RandomVariation]);
                     Variations = setfield(Variations,Fields{k},RandomVariation(1)*ones(1,nTissues));
                  end
               end
        end
        
        function Tp = calcTissueParameter(this,cst,numOfVoxels,ctScen)
               TissueParameters = this.BioModels{1}.calcTissueParameter(cst,numOfVoxels,ctScen);
               if this.includeNominalScenario
                  Tp(1) = {TissueParameters};

                  for k=2:this.nSamples
                     Tp(k) = {this.calcParameterVariations(TissueParameters)};
                  end
               else
                  for k=1:this.nSamples
                     Tp(k) = {this.calcParameterVariations(TissueParameters)};
                  end

               end

        end

        function set.nSamples(this,nSamples)
            valid = isnumeric(nSamples) && isscalar(nSamples) && mod(nSamples,1) == 0 && nSamples > 0;
            matRad_cfg = MatRad_Config.instance();
            if ~valid 
                matRad_cfg.dispError('Invalid value for nSamples! Needs to be a positive integer!');
            end
            this.nSamples = nSamples;
            this.updateScenarios();
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
            
            if ~isempty(this.SelectedBioModels)
               this.updateBioModel(this.SelectedBioModels);
            end
            
            
        end
        
        
        function this = updateBioModel(this,BioModel)
           BioModelInstance = this.setBioModel(BioModel{1});
           this.NumModelParameters = BioModelInstance{1}.NumModelParameters;
           this.BioModels = BioModelInstance;
           this.ModelParameters = BioModelInstance{1}.ModelParameters;
        end
        
        
    end


end