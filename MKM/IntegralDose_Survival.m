classdef IntegralDose_Survival < handle
   properties
      zD;
      survivalExecutionCommand;
      wDir;

      survivalSourcePath;
      
      survivalParameterFileName;
      calcProperties;
      wDirWsl;
      zD_singleIon;
   end
   
   methods
      function obj = IntegralDose_Survival()
         
      end
      
      function genParameterFile(obj)
         
         if ~exist(obj.wDir, 'dir')
            mkdir(obj.wDir);
         end
      
         fileName = strcat(obj.wDir, filesep, obj.survivalParameterFileName, '.txt');
         fID = fopen(fileName, 'w');
         %source
         fprintf(fID, 'source ');
         fprintf(fID, obj.survivalSourcePath);
         fprintf(fID, strcat('/setenv.sh \n'));       
         fprintf(fID, '\n');
         
         %projectName
         fprintf(fID, 'projectName="');
         fprintf(fID,obj.calcProperties.projectName);
         fprintf(fID, '"\n');
         fprintf(fID, '\n');
         
         %output
         fprintf(fID, 'output="');
         fprintf(fID,obj.calcProperties.output);
         fprintf(fID, '"\n');
         fprintf(fID, '\n');
         
         %model
         fprintf(fID, 'model="');
         fprintf(fID,obj.calcProperties.model);
         fprintf(fID, '"\n');
         
         %calculusType
         fprintf(fID, 'calculusType="');
         fprintf(fID,obj.calcProperties.calculusType);
         fprintf(fID, '"\n');
         fprintf(fID, '\n');
         
         %presicision
%          fprintf(fID, 'precision=%1.2f"',obj.calcProperties.precision );
%          fprintf(fID, '"\n');
%          fprintf(fID, '\n');
         
         %parallelismType
         fprintf(fID, 'parallelismType="');
         fprintf(fID,obj.calcProperties.parallelismType);
         fprintf(fID, '" \n');
         fprintf(fID, '\n');
         
         %cellType
         fprintf(fID, 'cellType="');
         fprintf(fID,obj.calcProperties.cellType);
         fprintf(fID, '" \n');
         fprintf(fID, '\n');
         
         %Model parameters
         switch obj.calcProperties.model
            case 'MKM'
 
               fprintf(fID, 'MKM_alpha0=%1.1f"', obj.calcProperties.modelParam.alpha0);
               fprintf(fID, '" \n');
               
               fprintf(fID, 'MKM_beta0=%1.1f"',obj.calcProperties.modelParam.beta0);
               fprintf(fID, '" \n');
               
               fprintf(fID, 'MKM_rNucleus=%1.1f"', obj.calcProperties.modelParam.rNucleus);
               fprintf(fID, '" \n');
               
               fprintf(fID, 'MKM_rDomain=%1.2f"',obj.calcProperties.modelParam.rDomain);

               fprintf(fID, '"\n');
               fprintf(fID, '\n');
         end
         
         %ion
         fprintf(fID, 'ion="');
         fprintf(fID,obj.calcProperties.ion);
         fprintf(fID, '"\n');
         fprintf(fID, '\n');
         
         %trackMode
         fprintf(fID, 'trackMode="');
         fprintf(fID,obj.calcProperties.trackMode);
         fprintf(fID, '"\n');
         fprintf(fID, '\n');
         
         %energies
         fprintf(fID, 'energies="');
         fprintf(fID, '%3.4f ', obj.calcProperties.energies);
         fprintf(fID, '"\n');
         fprintf(fID, '\n');
         
         %call
         fprintf(fID, ['cd ', obj.wDirWsl, '/\n']);
         fprintf(fID, '\n');
         fprintf(fID, strcat(obj.survivalSourcePath,  '/survival'));
         fields = fieldnames(obj.calcProperties);
         for k=1:size(fields,1)
            if ~strcmp(fields{k}, 'modelParam')
               fprintf(fID, [' -', fields{k}, ' $', fields{k}, ' \\\n']);
               fprintf(fID, '\t \t');
            end
         end

         fclose(fID);
      end
      
      
      function execute(obj)  
         idx = strfind(obj.survivalExecutionCommand, '\');
         executionCommand = obj.survivalExecutionCommand;

         executionCommand(idx) = '/';
         system(executionCommand, '-echo');
         
      end
      
      
      function readSingleIonLUT(obj,Ion)

         DataTable = importfileSurvival(strcat(obj.wDir, filesep, 'LUT_', num2str(Ion),'_LQparameters_MKM.csv'));
         obj.zD_singleIon.E = DataTable.meanEnergy;
         obj.zD_singleIon.zD = DataTable.alpha./DataTable.beta;

      end
      
      function plotZDSingleIon(obj,HOLD)
         if ~HOLD
            figure;
         else
            hold on;
         end
         semilogx(obj.zD_singleIon.E, obj.zD_singleIon.zD, '.-');
         grid on;
         hold off;
      end
      
   end
end