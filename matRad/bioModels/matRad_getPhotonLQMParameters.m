function [ax,bx,cellLines,cellLineIdx] = matRad_getPhotonLQMParameters(cst,numVoxel,VdoseGrid)
% matRad function to receive the photon LQM reference parameter
%
% call:
%   [ax,bx,cellLineNames,cellLineIdx] = matRad_getPhotonLQMParameters(cst,numVoxel,VdoseGrid)
%
% input:
%   cst:            matRad cst struct. cst{i,5}.bioParameters is expected
%                   to hold the fields cellLine (char), alphaX and betaX
%   numVoxel:       number of voxels of the dose cube
%   VdoseGrid:      optional linear index vector that allows to specify subindices
%                   for which ax and bx will be computed
%
% output:
%   ax:             vector containing for each linear voxel index alpha_x
%   bx:             vector containing for each linear voxel index beta_x
%   cellLineNames:  cell array of the unique cell line names found in cst
%   cellLineIdx:    vector containing for each linear voxel index the
%                   index into cellLineNames of the assigned cell line
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018-2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();

numOfCtScen = unique(cellfun(@numel,cst(:,4)));

if numel(numOfCtScen) > 1
    matRad_cfg.dispError('Inconsinstent number of ct scnearios in cst!');
end

for i = 1:size(cst,1)
    if isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET')

        if ~isfield(cst{i,5}, 'bioParameters') && isfield(cst{i,5}, 'alphaX') && isfield(cst{i,5}, 'betaX')
            matRad_cfg.dispDeprecationWarning('alphaX and betaX propoerties in cst will be deprecated soon! Use bioParameters instead.');
            % TODO: this needs to change to avoid cellLines with the same
            % name in dij.cellLineNames
            cst{i,5}.bioParameters.cellLine = 'Generic';
            cst{i,5}.bioParameters.alphaX   = cst{i,5}.alphaX;
            cst{i,5}.bioParameters.betaX    = cst{i,5}.betaX;
        end
    end
end

relevantIx = strcmp(cst(:,3),'OAR') | strcmp(cst(:,3),'TARGET');

cellLineNames = cellfun(@(voi) char(voi.bioParameters.cellLine), cst(relevantIx,5), 'UniformOutput', false);
[cellLineNames, structIdx] = unique(cellLineNames);

cellLineAx = cellfun(@(voi) voi.bioParameters.alphaX, cst(relevantIx,5), 'UniformOutput', false);
cellLineAx = cellLineAx(structIdx);

cellLineBx = cellfun(@(voi) voi.bioParameters.betaX, cst(relevantIx,5), 'UniformOutput', false);
cellLineBx = cellLineBx(structIdx);

ax = cell(numOfCtScen,1);
bx = cell(numOfCtScen,1);
cellLineIdx = cell(numOfCtScen,1);

[ax{:}] = deal(zeros(numVoxel,1));
[bx{:}] = deal(zeros(numVoxel,1));
[cellLineIdx{:}] = deal(zeros(numVoxel,1));

for i = 1:size(cst,1)
    if isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET')

        bioParameters = cst{i,5}.bioParameters;
        cellLineIx = find(strcmp(cellLineNames,bioParameters.cellLine));

        for s = 1:numOfCtScen

            if exist('VdoseGrid','var')
                if iscell(VdoseGrid)
                    isInVdoseGrid = ismember(VdoseGrid{s},cst{i,4}{s});
                else
                    isInVdoseGrid = ismember(VdoseGrid,cst{i,4}{s});
                end
                ax{s}(VdoseGrid(isInVdoseGrid)) = bioParameters.alphaX;
                bx{s}(VdoseGrid(isInVdoseGrid)) = bioParameters.betaX;
                cellLineIdx{s}(VdoseGrid(isInVdoseGrid)) = cellLineIx;


            else
                ax{s}(cst{i,4}{s}) = bioParameters.alphaX;
                bx{s}(cst{i,4}{s}) = bioParameters.betaX;
                cellLineIdx{s}(cst{i,4}{s}) = cellLineIx;
            end

        end
    end
end


for cLine=1:numel(cellLineNames)

    cellLines(cLine).cellLineName = cellLineNames{cLine};
    cellLines(cLine).alphaX       = cellLineAx{cLine};
    cellLines(cLine).betaX       = cellLineBx{cLine};
    
end