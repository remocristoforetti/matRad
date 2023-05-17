function model = matRad_bioModel(sRadiationMode,sQuantityOpt, sModel, machine)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  matRad_biosModel
%  This is a helper function to instantiate a matRad_BiologicalsModel. This
%  function currently exists for downwards compatability, as the new
%  Biological sModels will follow a polymorphic software architecture
%
% call
%   matRad_biosModel(ssRadiationMode,sQuantityOpt, ssModel)
%
%   e.g. pln.bioParam = matRad_biosModel('protons','constRBE','RBExD')
%
% input
%   ssRadiationMode:     radiation modality 'photons' 'protons' 'carbon'
%   sQuantityOpt:       string to denote the quantity used for
%                       optimization 'physicalDose', 'RBExD', 'effect'
%   ssModel:             string to denote which biological sModel is used
%                       'none': for photons, protons, carbon                'constRBE': constant RBE for photons and protons
%                       'MCN': McNamara-variable RBE sModel for protons      'WED': Wedenberg-variable RBE sModel for protons
%                       'LEM': Local Effect sModel for carbon ions
%
% output
%   sModel:              instance of a biological sModel
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%sModel = matRad_BiologicalsModel(sRadiationMode,sQuantityOpt,sModel);
    switch sModel
         case{'none'}
             if (strcmp(sRadiationMode, 'photons'))
                 model = matRad_bioModel_constRBE(sRadiationMode);
             else
                model = matRad_bioModel_noneParticles(sRadiationMode);
             end

         case{'constRBE'}
             model = matRad_BioModel_constRBE(sRadiationMode);
%         %For protons only, might be set as a further subclass, all require
%         %LET
        case {'MCN'}
             model = matRad_bioModelMCN(sRadiationMode);
        case {'WED'}
             model = matRad_bioModelWED(sRadiationMode);
        case {'ESTRO'}
             model = matRad_bioModelESTRO(sRadiationMode);
        case {'LSM'}
             model = matRad_BioModel_LSM(sRadiationMode);
        case {'MKMLET'}
             model = matRad_bioModelMKMLET(sRadiationMode);
        case {'MKMLET_corrected'}
             model = matRad_BioModel_MKMLET_corrected(sRadiationMode);
        case{'CAR'}
             model = matRad_bioModelCAR(sRadiationMode);
%         %for heavier particles

       case {'LEM'}
            model = matRad_BioModel_LEM(sRadiationMode);
        case {'MKM_Kase2008'}
            model = matRad_BioModel_SurvivalMKM_Kase2008(sRadiationMode);
       case {'tabulatedRBEModel'}
            model = matRad_tabulatedRBEModel(sRadiationMode);
       otherwise
           matRad_cfg =  MatRad_Config.instance();
           matRad_cfg.dispWarning('This sModel has not yet been implemented \n');
           %Set empty generic sModel ?
           if strcmp(sRadiationMode, 'photons')
                model = matRad_BioModel_constRBE(sRadiationMode);
           else
                model = matRad_BioModel_noneParticles(sRadiationMode);   
           end
    end
    
    if ~isempty(model)
        model.quantityOpt = sQuantityOpt;
        if nargin>3 && ~isempty(model.RequiredBaseData)
         model.baseData = machine;
         
        end
    end

end % end class definition