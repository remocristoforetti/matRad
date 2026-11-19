%% Script: Write a TOPAS Schneider converter with the MCsquare (matRad_default) HU-to-density relation
% TOPAS computes the voxel density as
%   rho = (Offset + Factor*(FactorOffset + HU)) * DensityCorrection(HU - MinImagingValue + 1)
% with one DensityCorrection entry per integer HU. Using a single section
% with Offset = 1 and Factor = 0 makes rho = DensityCorrection(HU), i.e. an
% explicit per-HU value, linearly interpolated from the MCsquare table
% (the same interpolation MCsquare applies internally). A single water
% material is assigned to all HU, as in MCsquare's matRad_default scanner.
clear;
matRad_rc;
matRad_cfg = MatRad_Config.instance();

mcSquareTable = fullfile(matRad_cfg.matRadRoot,'thirdParty','MCsquare','bin','Scanners','matRad_default','HU_Density_Conversion.txt');
outFile       = fullfile(matRad_cfg.primaryUserFolder,'dataPaper','doseCalculation','MCengines','TOPAS_SchneiderConverter_MCsquareDensity.txt');

txt = regexprep(fileread(mcSquareTable),'#[^\n]*','');
tab = reshape(sscanf(txt,'%f'),2,[])';

huMin = -1024;
huMax = 3071;
hu    = (huMin:huMax)';
rho   = interp1(tab(:,1),tab(:,2),hu,'linear','extrap');

fid = fopen(outFile,'w');
fprintf(fid,'# Schneider converter with the MCsquare matRad_default HU-to-density relation (one density per HU)\n');
fprintf(fid,'# generated from %s\n\n',strrep(mcSquareTable,'\','/'));
fprintf(fid,'i:Ge/Patient/MinImagingValue = %d\n',huMin);
fprintf(fid,['dv:Ge/Patient/DensityCorrection = %d' repmat(' %.6g',1,numel(rho)) ' g/cm3\n'],numel(rho),rho);
fprintf(fid,'iv:Ge/Patient/SchneiderHounsfieldUnitSections = 2 %d %d\n',huMin,huMax + 1);
fprintf(fid,'uv:Ge/Patient/SchneiderDensityOffset = 1 1.0\n');
fprintf(fid,'uv:Ge/Patient/SchneiderDensityFactor = 1 0.0\n');
fprintf(fid,'uv:Ge/Patient/SchneiderDensityFactorOffset = 1 0.0\n\n');
fprintf(fid,'sv:Ge/Patient/SchneiderElements = 2 "Hydrogen" "Oxygen"\n');
fprintf(fid,'iv:Ge/Patient/SchneiderHUToMaterialSections = 2 %d %d\n',huMin,huMax + 1);
fprintf(fid,'uv:Ge/Patient/SchneiderMaterialsWeight1 = 2 0.111894 0.888106\n');
fprintf(fid,'dv:Ge/Patient/SchneiderMaterialMeanExcitationEnergy = 1 78.0 eV\n');
fclose(fid);

matRad_cfg.dispInfo('TOPAS Schneider converter written to %s\n',outFile);
