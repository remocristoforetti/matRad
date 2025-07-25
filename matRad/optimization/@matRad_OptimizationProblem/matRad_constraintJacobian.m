function jacob = matRad_constraintJacobian(optiProb,w,dij,cst)
% matRad IPOPT callback: jacobian function for inverse planning 
% supporting max dose constraint, min dose constraint, min mean dose constraint, 
% max mean dose constraint, min EUD constraint, max EUD constraint, max DVH 
% constraint, min DVH constraint 
% 
% call
%   jacob = matRad_jacobFunc(optiProb,w,dij,cst)
%
% input
%   optiProb: option struct defining the type of optimization
%   w:    bixel weight vector
%   dij:  dose influence matrix
%   cst:  matRad cst struct
%
% output
%   jacob: jacobian of constraint function
%
% References
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get current dose / effect / RBExDose vector
%d = matRad_backProjection(w,dij,optiProb);
%d = optiProb.matRad_backProjection(w,dij);
optiProb.BP.compute(dij,w);
d = optiProb.BP.GetResult();

% initialize jacobian (only single scenario supported in optimization)
% jacob = sparse([]);
% 
% % initialize projection matrices and id containers
% DoseProjection{1}          = sparse([]);
% mAlphaDoseProjection{1}    = sparse([]);
% mSqrtBetaDoseProjection{1} = sparse([]);
% voxelID                     = [];
% constraintID                = [];
% 
% % get the used scenarios
% useScen  = optiProb.BP.scenarios;
% scenProb = optiProb.BP.scenarioProb;

% retrieve matching 4D scenarios
% fullScen      = cell(ndims(d),1);
% [fullScen{:}] = ind2sub(size(d),useScen);
% contourScen   = fullScen{1};

% compute objective function for every VOI.
constIndexigStruct = [optiProb.constrIdx, zeros(size(optiProb.constrIdx,1),1)];
for i = 1:size(optiProb.constrIdx,1)
   
      constraint = optiProb.constraints{i}; %Get the Optimization Object
      curConIdx = constIndexigStruct(i,1); %optiProb.constrIdx(i,1);

      robustness = constraint.robustness;

      if isa(constraint, 'DoseConstraints.matRad_DoseConstraint')
        
          quantityConstrained = constraint.quantity;
          quantityNames = cellfun(@(x) x.quantityName,optiProb.BP.quantities, 'UniformOutput',false);
          quantityConstrainedInstance = optiProb.BP.quantities{strcmp(quantityConstrained,quantityNames)};
          
            if ~exist('fJacob', 'var') || ~isfield(fJacob,quantityConstrained)
                if isa(quantityConstrainedInstance, 'matRad_DistributionQuantity')
                    fJacob.(quantityConstrained) = cell(1);
                    %fJacob.(quantityConstrained) = {zeros(dij.doseGrid.numOfVoxels,1)};
                elseif isa(quantityConstrainedInstance, 'matRad_ScalarQuantity')
                    fJacob.(quantityConstrained)= cell(size(cst,1),1);
                end
            
            end
          
          switch robustness
             
             case 'none' % if conventional opt: just sum objectiveectives of nominal dose
                  if isa(quantityConstrainedInstance, 'matRad_DistributionQuantity')
                      structIdxs = cst{curConIdx,4}{1};
                      d_i = d.(quantityConstrained){1}(structIdxs);
                  elseif isa(quantityConstrainedInstance, 'matRad_ScalarQuantity')
                      d_i = d.(quantityConstrained){curConIdx};
                  end

                  jacobSub = constraint.computeDoseConstraintJacobian(d_i);
                      
             case 'PROB' % if prob opt: sum up expectation value of objectives
                  if isa(quantityConstrainedInstance, 'matRad_DistributionQuantity')
                      structIdxs = cat(1,cst{curConIdx,4}{:});
                      structIdxs = unique(structIdxs);

                      d_i = d.(quantityConstrained){1}(structIdxs);
                  elseif isa(quantityConstrainedInstance, 'matRad_ScalarQuantity')
                      d_i = d.(quantityConstrained){curConIdx};
                  end
                
                jacobSub = constraint.computeDoseConstraintJacobian(d_i);
                
             case 'VWWC'  % voxel-wise worst case - takes minimum dose in TARGET and maximum in OAR
                contourIx = unique(contourScen);
                if ~isscalar(contourIx)
                      % voxels need to be tracked through the 4D CT,
                      % not yet implemented
                      matRad_cfg.dispError('4D VWWC optimization is currently not supported');
                end
                
                % prepare min/max dose vector
                if ~exist('d_tmp','var')
                      d_tmp = [d{useScen}];
                end
                
                d_Scen = d_tmp(cst{curConIdx,4}{contourIx},:);
                
                d_max = max(d_Scen,[],2);
                d_min = min(d_Scen,[],2);
                
                if isequal(cst{curConIdx,3},'OAR')
                      d_i = d_max;
                elseif isequal(cst{curConIdx,3},'TARGET')
                      d_i = d_min;
                end
                jacobSub = constraint.computeDoseConstraintJacobian(d_i);
                
             case 'VWWC_INV'  %inverse voxel-wise conformitiy - takes maximum dose in TARGET and minimum in OAR
                contourIx = unique(contourScen);
                if ~isscalar(contourIx)
                      % voxels need to be tracked through the 4D CT,
                      % not yet implemented
                      matRad_cfg.dispError('4D inverted VWWC optimization is currently not supported');
                end
                
                % prepare min/max dose vector
                if ~exist('d_tmp','var')
                      d_tmp = [d{useScen}];
                end
                
                d_Scen = d_tmp(cst{curConIdx,4}{contourIx},:);
                
                d_max = max(d_Scen,[],2);
                d_min = min(d_Scen,[],2);
                
                if isequal(cst{curConIdx,3},'OAR')
                      d_i = d_min;
                elseif isequal(cst{curConIdx,3},'TARGET')
                      d_i = d_max;
                end
                jacobSub = constraint.computeDoseConstraintJacobian(d_i);
                
             otherwise
                matRad_cfg.dispError('Robustness setting %s not yet supported!',constraint.robustness);
          end
          
          nConst = size(jacobSub,2);
          constIndexigStruct(i,3) = nConst;
            if isa(quantityConstrainedInstance, 'matRad_DistributionQuantity')
                
                startIx = size(fJacob.(quantityConstrained){1},2) + 1;
                fJacob.(quantityConstrained) = {[fJacob.(quantityConstrained){1},sparse(dij.doseGrid.numOfVoxels,nConst)]};
                fJacob.(quantityConstrained){1}(structIdxs,startIx:end) = jacobSub;
            
            elseif isa(quantityConstrainedInstance, 'matRad_ScalarQuantity')
                startIx = size(fJacob.(quantityConstrained){curConIdx},2) + 1;

                fJacob.(quantityConstrained)(curConIdx) = {[fJacob.(quantityConstrained){curConIdx},sparse(1,nConst)]}; 
                fJacob.(quantityConstrained){curConIdx}(1,startIx:end) = jacobSub;
            end

      elseif isa(constraint, 'OmegaConstraints.matRad_VarianceConstraint')

         quantityConstrained = constraint.quantity;
         quantityNames = cellfun(@(x) x.quantityName,optiProb.BP.quantities, 'UniformOutput',false);
         quantityConstrainedInstance = optiProb.BP.quantities{strcmp(quantityConstrained,quantityNames)};

         if ~exist('fJacob', 'var') || ~isfield(fJacob,quantityConstrained)
            if isa(quantityConstrainedInstance, 'matRad_DistributionQuantity')
                fJacob.(quantityConstrained) = cell(1);
            elseif isa(quantityConstrainedInstance, 'matRad_ScalarQuantity')
                fJacob.(quantityConstrained) = cell(size(cst,1),1);           
            end
        
         end
         
         if isa(quantityConstrainedInstance, 'matRad_CumulativeScalarQuantity')
            d_i = d.(quantityConstrained){curConIdx};
            jacobSub = constraint.computeVarianceConstraintJacobian(d_i, 1);    
         else
             allVoxels = arrayfun(@(scenStruct) scenStruct{1}, cst{curConIdx,4}, 'UniformOutput',false);
             nVoxels = numel(unique(vertcat(allVoxels{:}))); 
    
    
             switch robustness
    
                case 'PROB'
                    d_i = d.(quantityConstrained){curConIdx};
                    jacobSub = constraint.computeVarianceConstraintJacobian(d_i, nVoxels);    
             end
         end

        nConst = size(jacobSub,2);
        
        constIndexigStruct(i,3) = nConst;
        
        %startIx = size(fJacob.(quantityConstrained){curConIdx},2) + 1;
        fJacob.(quantityConstrained){curConIdx} = [fJacob.(quantityConstrained){curConIdx},sparse(jacobSub)];

      end
end

for qtIdx=optiProb.BP.constrainedQuantities
    % this might cause problems with distribution quantities
    fJacob.(qtIdx{1})(find(cellfun(@(x) isempty(x), fJacob.(qtIdx{1})))) = {zeros(size(find(cellfun(@(x) ~isempty(x), fJacob.(qtIdx{1})), 1, 'first'),1))};
end

optiProb.BP.computeConstraintJacobian(dij,fJacob,w);

j = optiProb.BP.wJacob;

jacob = [];

cumConstrDistIdx = struct();
cumConstrScalarIdx = cell(size(cst,1),1);
for i = 1:size(optiProb.constrIdx,1)
    constraint = optiProb.constraints{i};
    
    quantityNames = cellfun(@(x) x.quantityName,optiProb.BP.quantities, 'UniformOutput',false);
    quantityConstrainedInstance = optiProb.BP.quantities{strcmp(constraint.quantity,quantityNames)};
    
    curStructIdx        = optiProb.constrIdx(i,1);
    %curConstInStructIdx = optiProb.constrIdx(i,2);
    if isa(quantityConstrainedInstance, 'matRad_DistributionQuantity')
        % Only one scenario for now. Need to pick the contraints entry in
        % the right order

        if ~isfield(cumConstrDistIdx, constraint.quantity)
            cumConstrDistIdx.(constraint.quantity) = 1;
        end
         
        currConstrDistIdx = cumConstrDistIdx.(constraint.quantity);
        % Collect the constraint jacobian for this quantity in the right
        % order
        jacob = [jacob; j.(constraint.quantity){1}(currConstrDistIdx:currConstrDistIdx + constIndexigStruct(i,3) -1,:)];

        cumConstrDistIdx.(constraint.quantity) = cumConstrDistIdx.(constraint.quantity) + constIndexigStruct(i,3);
    elseif isa(quantityConstrainedInstance, 'matRad_ScalarQuantity')
        
        if ~isfield(cumConstrScalarIdx{curStructIdx}, constraint.quantity)
            cumConstrScalarIdx{curStructIdx}.(constraint.quantity) = 1;
        end
        
        currConstrScalarIdx = cumConstrScalarIdx{curStructIdx}.(constraint.quantity);

        jacob = [jacob; j.(constraint.quantity){curStructIdx}(currConstrScalarIdx:currConstrScalarIdx + constIndexigStruct(i,3) -1,:)];

        cumConstrScalarIdx{curStructIdx}.(constraint.quantity) = cumConstrScalarIdx{curStructIdx}.(constraint.quantity) + constIndexigStruct(i,3);

    end
end
% constrQts = optiProb.BP.constrainedQuantities;
% [~,orderedFields] = intersect(fieldnames(fJacob)', constrQts);
% 
% jacob = sparse([]);
% for qtIdx=fieldnames(j)'%optiProb.BP.constrainedQuantities
%     nScensOrStructs   = find(cellfun(@(x) ~isempty(x) && nnz(x) ~= 0, j.(qtIdx{1})))';
%     for elementIdx=nScensOrStructs
%         jacob = [jacob; j.(qtIdx{1}){elementIdx}];
%     end
% end

    
gradientChecker = 0;
if gradientChecker == 1
    f =  matRad_constraintFunctions(optiProb,w,dij,cst);
    epsilon = 1e-5;

    ix = randi([dij.totalNumOfBixels],numel(f),5);

    for jIx=1:numel(f)
        for i=ix(jIx,:)
    
            wInit = w;
            wInit(i) = wInit(i) + epsilon;
            fDel = matRad_constraintFunctions(optiProb,wInit,dij,cst);
            numGrad = (fDel(jIx) - f(jIx))/epsilon;
            diff = (numGrad/jacob(jIx,i) - 1)*100;
            fprintf(['grad val #' num2str(i) '- rel diff numerical and analytical jacobian = ' num2str(diff) '\n']);
            %fprintf([' any nan or zero for photons' num2str(sum(isnan(glog{1}))) ',' num2str(sum(~logical(glog{1}))) ' for protons: ' num2str(sum(isnan(glog{2}))) ',' num2str(sum(~logical(glog{2}))) '\n']);
        end
    end
end


end
