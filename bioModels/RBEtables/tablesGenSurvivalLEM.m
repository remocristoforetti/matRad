%Generate RBEtables interfacing with survival code
matRad_cfg = MatRad_Config.instance();
%Define energies (MeV/u)
energies = logspace(-1,log10(300),250);


%Define ions
ions  = {'H', 'He', 'Li', 'Be', 'B', 'C'};
ionsMassNumber = [1 4 7 9 11 12]; 

%Define model and calculus type
model        = 'LEMI';
calculusType = 'rapidLEM_Russo2011';

%Define parameters to compute (for the time being fix alpha/beta)

alpha0 = [0.313];

beta0  = [0.062];

%MKM specific
Dt     = [19];                      %Taken from LEM3 paper, figure 5 (ACCURACY OF THE LOCAL EFFECT MODEL FOR THE PREDICTION OF BIOLOGIC EFFECTS OF CARBON ION BEAMS IN VITRO AND IN VIVO)
RN     = 5*ones(1,size(Dt,2));

%Define the interface class
IntegralS = integralDoseSurvival();
IntegralS.wDir = [matRad_cfg.matRadRoot, filesep, 'survival'];

%This is source path see from wsl
IntegralS.survivalSourcePath = '/mnt/c/Users/r408i/Desktop/r408i_data/Survival';

IntegralS.survivalParameterFileName = 'LEM3';
IntegralS.calcProperties.output          = 'LQ_pars';
IntegralS.calcProperties.model           = model;
IntegralS.calcProperties.calculusType    = calculusType;
%IntegralS.calcProperties.precision       = 0.15;
IntegralS.calcProperties.parallelismType = '0';
IntegralS.calcProperties.cellType        = 'HSG';

IntegralS.calcProperties.trackMode = 'histogram';
IntegralS.calcProperties.energies  = energies;

idx = strfind(IntegralS.wDir, '\');
dirPath = IntegralS.wDir(idx(1)+1:end);
dirPath(strfind(dirPath, '\')) = '/';
IntegralS.wDirWsl = ['/mnt/c/',dirPath];


for alphaIdx = 1:size(alpha0,2)
      IntegralS.calcProperties.modelParam.LEM_alpha0   = alpha0(alphaIdx);
      IntegralS.calcProperties.modelParam.LEM_beta0    = beta0(alphaIdx);
    
      for pIdx=1:size(Dt,2)
        
        IntegralS.calcProperties.projectName     = ['LUT_alpha_',num2str(alpha0(alphaIdx)),'_Dt_', num2str(Dt(pIdx))];
      
        IntegralS.calcProperties.modelParam.LEM_rNucleus = RN(pIdx);
        IntegralS.calcProperties.modelParam.LEM_Dt  = Dt(pIdx);
        matRad_cfg.dispInfo('alpha value %d, Dt value %d \n', alpha0(alphaIdx), Dt(pIdx));
        for k=1:size(ions,2)
            IntegralS.calcProperties.energies  = energies.*ionsMassNumber(k);
            IntegralS.calcProperties.ion = ions{k};
            IntegralS.genParameterFile();
            IntegralS.survivalExecutionCommand = ['wsl ',IntegralS.wDirWsl,'/', IntegralS.survivalParameterFileName, '.txt'];
        
            %Execution with system does not work, check -> Solved, if does not work
            %again, check that default destribution is Ubuntu and not Ubunutu-20.04.
            %Else, from prompt wsl --setdefault Ubuntu
            matRad_cfg.dispInfo('Executing Survival...');
            a = IntegralS.execute();
            if a == 0
                matRad_cfg.dispInfo('done. \n');
            else
                matRad_cfg.dispInfo('error. \n');
            end
        end
      end

end

%% read Out
clear RBEtable;
%meta info
for pIdx = 1:size(Dt,2)
    RBEtable.meta.model = [calculusType, '_Dt_', num2str(Dt(pIdx))];
    RBEtable.meta.description = ['Model obtained with', calculusType, ' computed with Survival code.'];
   
    modelParameters = struct('alphaX', [], ...
                             'betaX', [], ...
                             'rNucleus', [], ...
                             'Dt', []);
                          
    for alphaIdx =1:size(alpha0,2)        
        modelParameters.alphaX = alpha0;
        modelParameters.betaX  = beta0;
        modelParameters.rNucleus = RN(pIdx);
        modelParameters.Dt  = Dt(pIdx);
    
        RBEtable.meta.modelParameters = modelParameters;


        %filenames = dir([IntegralS.wDir, filesep, '*.csv']);
        filenames = [IntegralS.wDir, filesep,'LUT_alpha_',num2str(alpha0(alphaIdx)),'_Dt_', num2str(Dt(pIdx)), '_LQparameters_LEMI.csv'];
        [alphaE, betaE] = IntegralS.readMultipleIonLUT(ions, filenames);
        
        %This is directly alpha, not zD;
%         for iIdx = 1:size(ions)
%             varAlpha(:,k) = alpha0(alphaIdx) + beta0(aphaIdx) * alphaE;
%             varBeta(:,k) = betaE;
%         end
        
        RBEtable.data(alphaIdx).alphaX = alpha0(alphaIdx);
        RBEtable.data(alphaIdx).betaX = beta0(alphaIdx);

        RBEtable.data(alphaIdx).energies     = energies;
        RBEtable.data(alphaIdx).includedIonZ = [1 2 3 4 5 6];
        for m=1:size(alphaE,2)

            RBEtable.data(alphaIdx).alpha(:,m)        = alphaE{m};
            RBEtable.data(alphaIdx).beta(:,m)         = betaE{m};
        end
%         RBEtable.data(alphaIdx).alpha        = alphaE;
%         RBEtable.data(alphaIdx).beta         = betaE;
    end

    save([matRad_cfg.matRadRoot, filesep, 'bioModels', filesep, 'RBEtables', filesep,'RBEtable_', calculusType, '_Paper_HSG_LEMI', num2str(Dt(pIdx)), '.mat'], 'RBEtable');
    RBEtable = [];
end