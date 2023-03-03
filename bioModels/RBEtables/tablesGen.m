%table data generator

%This is just an example, load existing lUT
RBEtable.meta.model = 'MKM_Inaniwa2010';
RBEtable.meta.description = 'Model obtained with Inaniwa 2010 correction for Zsat, computed with Matlab code. Example';


modelParameters = struct('alphaX', [], ...
                         'betaX', [], ...
                         'rNucleus', [], ...
                         'rDomain', []);
                      
                      
modelParameters.alphaX = [0.05, 0.1, 0.15];
modelParameters.betaX  = [0.05, 0.05, 0.05];
modelParameters.rNucleus = 3.9;
modelParameters.rDomain  = 0.35;

RBEtable.meta.modelParameters = modelParameters;

%table.data.energies = linspace(0.1,300,100)';

pathL = 'D:\matRad_gitHubRemo_MKM_TOPAS\MKM\LUTs_rn35_RN_39_Beta005_500MeV_opt';
filenames = dir([pathL, filesep,'*.mat']);
for k=1:size(filenames,1)

   ZD(k) = load([filenames(1).folder, filesep, filenames(k).name]);
   
end

for tissueIdx = 1:size(modelParameters.alphaX,2)
   for k=1:size(ZD,2)
      alpha(:,k) = modelParameters.alphaX(tissueIdx) + (ZD(k).zD).*modelParameters.betaX(tissueIdx);
      beta(:,k)  = modelParameters.betaX(tissueIdx)*ones(size(ZD(k).E,2),1);
   end
   
   RBEtable.data(tissueIdx).alphaX = modelParameters.alphaX(tissueIdx);
   RBEtable.data(tissueIdx).betaX = modelParameters.betaX(tissueIdx);
   
   RBEtable.data(tissueIdx).energies = ZD(1).E';
   RBEtable.data(tissueIdx).includedIonZ = [1 2 3 4 5 6];
   RBEtable.data(tissueIdx).alpha = alpha;
   RBEtable.data(tissueIdx).beta  = beta;

end

%table.data.alpha = 

