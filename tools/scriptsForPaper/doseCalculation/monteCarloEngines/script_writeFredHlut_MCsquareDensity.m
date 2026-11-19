%% Script: Write a FRED HU table with the MCsquare (matRad_default) HU-to-density relation
% FRED interpolates linearly between 'mat:' rows, so the MCsquare
% breakpoints are copied one-to-one. All rows are water (as in MCsquare's
% matRad_default scanner, I = 78 eV), so rho = RSP. The table is written
% to userdata/hluts, where the FRED engine finds it by name:
%   engineFred.HUtable = 'matRad_MCsquareDensity_FredMaterialConverter';
clear;
matRad_rc;
matRad_cfg = MatRad_Config.instance();

mcSquareTable = fullfile(matRad_cfg.matRadRoot,'thirdParty','MCsquare','bin','Scanners','matRad_default','HU_Density_Conversion.txt');
outDir        = fullfile(matRad_cfg.primaryUserFolder,'hluts');
outFile       = fullfile(outDir,'matRad_MCsquareDensity_FredMaterialConverter.txt');

txt = regexprep(fileread(mcSquareTable),'#[^\n]*','');
tab = reshape(sscanf(txt,'%f'),2,[])';

if ~isfolder(outDir)
    mkdir(outDir);
end
fid = fopen(outFile,'w');
fprintf(fid,'matColumns: HU rho RSP Ipot Lrad C Ca H N O P Ti S\n');
fprintf(fid,'mat: %g %.6g %.6g 78.0 36.1 0 0 11.189400 0 88.810600 0 0 0\n',[tab(:,1) tab(:,2) tab(:,2)]');
fclose(fid);

matRad_cfg.dispInfo('FRED HU table written to %s\n',outFile);
