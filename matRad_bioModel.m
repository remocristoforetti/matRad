function [bioModel] = matRad_bioModel(radiationMode,quantityOpt,model,varargin)

    switch model
         case{'none'}
             if (strcmp(radiationMode, 'photons'))
                 bioModel = matRad_BioModel_constRBE(radiationMode);
             else
                bioModel = matRad_BioModel_noneParticles(radiationMode);
             end

         case{'constRBE'}
             bioModel = matRad_BioModel_constRBE(radiationMode);
%         %For protons only, might be set as a further subclass, all require
%         %LET
        case {'MCN'}
             bioModel = matRad_BioModel_MCN(radiationMode);
        case {'WED'}
             bioModel = matRad_BioModel_WED(radiationMode);
        case {'LSM'}
             bioModel = matRad_BioModel_LSM(radiationMode);
        case {'MKMLET'}
             bioModel = matRad_BioModel_MKMLET(radiationMode);
        case {'MKMLET_corrected'}
             bioModel = matRad_BioModel_MKMLET_corrected(radiationMode);
        case{'CAR'}
             bioModel = matRad_BioModel_CAR(radiationMode);
%         %for heavier particles
        case {'LEM'}
            bioModel = matRad_BioModel_LEM(radiationMode);
       otherwise
           matRad_cfg =  MatRad_Config.instance();
           matRad_cfg.dispWarning('This model has not yet been implemented \n');
           %Set empty generic model ?
           if strcmp(radiationMode, 'photons')
                bioModel = matRad_BioModel_constRBE(radiationMode);
           else
                bioModel = matRad_BioModel_noneParticles(radiationMode);   
           end
    end
    
    if ~isempty(bioModel)
        bioModel.quantityOpt = quantityOpt;
        if nargin>3 && ~isempty(bioModel.RequiredBaseData)
         bioModel.baseData = varargin{:};
        end
    end
end