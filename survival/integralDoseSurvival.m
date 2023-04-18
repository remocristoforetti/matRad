classdef integralDoseSurvival < handle
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
      function obj = integralDoseSurvival()
         
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
               fprintf(fID, 'MKM_alpha0=%1.2f', obj.calcProperties.modelParam.MKM_alpha0);
               fprintf(fID, '\n');
               
               fprintf(fID, 'MKM_beta0=%1.2f',obj.calcProperties.modelParam.MKM_beta0);
               fprintf(fID, '\n');
               
               fprintf(fID, 'MKM_rNucleus=%1.2f', obj.calcProperties.modelParam.MKM_rNucleus);
               fprintf(fID, '\n');
               
               fprintf(fID, 'MKM_rDomain=%1.2f',obj.calcProperties.modelParam.MKM_rDomain);
               fprintf(fID, '\n');
               fprintf(fID, '\n');
             
             case {'LEMI', 'LEMII', 'LEMIII'}
               fprintf(fID, 'LEM_alpha0=%1.2f', obj.calcProperties.modelParam.LEM_alpha0);
               fprintf(fID, '\n');
               
               fprintf(fID, 'LEM_beta0=%1.2f',obj.calcProperties.modelParam.LEM_beta0);
               fprintf(fID, '\n');
               
               fprintf(fID, 'LEM_rNucleus=%1.2f', obj.calcProperties.modelParam.LEM_rNucleus);
               fprintf(fID, '\n');
               
               fprintf(fID, 'LEM_Dt=%1.2f',obj.calcProperties.modelParam.LEM_Dt);
               fprintf(fID, '\n');
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
         modelParametersNames = fieldnames(obj.calcProperties.modelParam);

         for k=1:size(fields,1)

             if ~strcmp(fields{k}, 'modelParam')
                fprintf(fID, [' -', fields{k}, ' $', fields{k}, ' \\\n']);
               fprintf(fID, '\t \t');
            else
               for m=1:size(modelParametersNames,1)
                fprintf(fID, [' -', modelParametersNames{m}, ' $', modelParametersNames{m}, ' \\\n']);
                fprintf(fID, '\t \t');

               end
            end
         end

         fclose(fID);
      end
      
      
      function a = execute(obj)
         idx = strfind(obj.survivalExecutionCommand, '\');
         executionCommand = obj.survivalExecutionCommand;

         executionCommand(idx) = '/';
         [a, ~] = system(executionCommand);
         
      end
      
      
      function [alphaE, betaE] = readSingleIonLUT(obj)

         DataTable = obj.importfileSurvival(strcat(obj.wDir, filesep, obj.calcProperties.projectName,'_LQparameters_MKM.csv'));
         obj.zD_singleIon.E = DataTable.meanEnergy;
         alphaE = DataTable.alpha;
         betaE = DataTable.beta;
         % Only valid if alpha0 = 0
         obj.zD_singleIon.zD = DataTable.alpha./DataTable.beta;

      end
      
      function [alphaE, betaE] = readMultipleIonLUT(obj,ions,filename)

         %DataTable = obj.importfileSurvival(strcat(obj.wDir, filesep, obj.calcProperties.projectName,'_LQparameters_MKM.csv'));
         %obj.zD_singleIon.E = DataTable.meanEnergy;
         DataTable = obj.importfileSurvival(filename);
         for k=1:size(ions,2)
            alphaE{k} = DataTable.alpha(DataTable.particle == ions{k});
            betaE{k} = DataTable.beta(DataTable.particle == ions{k});
         end
       end

%        function plotSingleLUT(obj,HOLD,alphaE)
%             if ~HOLD
%                 figure;
%              else
%                 hold on;
%             end
%             semilogx
%        end
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
      
      function KMKLUTLQparametersMKM = importfileSurvival(obj,filename, dataLines)
      %IMPORTFILE Import data from a text file
      %  KMKLUTLQPARAMETERSMKM = IMPORTFILE(FILENAME) reads data from text
      %  file FILENAME for the default selection.  Returns the data as a table.
      %
      %  KMKLUTLQPARAMETERSMKM = IMPORTFILE(FILE, DATALINES) reads data for
      %  the specified row interval(s) of text file FILENAME. Specify
      %  DATALINES as a positive scalar integer or a N-by-2 array of positive
      %  scalar integers for dis-contiguous row intervals.
      %
      %  Example:
      %  KMKLUTLQparametersMKM = importfile("C:\Users\Remo\Desktop\KMKLUT_LQparameters_MKM.csv", [2, Inf]);
      %
      %  See also READTABLE.
      %
      % Auto-generated by MATLAB on 17-Feb-2023 10:33:57

      %% Input handling

      % If dataLines is not specified, define defaults
      if nargin < 3
         dataLines = [2, Inf];
      end

      %% Set up the Import Options and import the data
      opts = delimitedTextImportOptions("NumVariables", 22);

      % Specify range and delimiter
      opts.DataLines = dataLines;
      opts.Delimiter = ",";

      % Specify column names and types
      opts.VariableNames = ["model", "calculusType", "cell", "alpha_0", "beta_0", "r_nucleus", "r_domain", "radiation_type", "particle", "meanEnergy", "meanEnergy_sigma", "meanLET", "meanLET_sigma", "LETd", "LETd_sigma", "nFractions", "timeSpacing", "fracDeliveryTime", "alpha", "beta", "alpha_sigma", "beta_sigma"];
      opts.VariableTypes = ["categorical", "double", "categorical", "double", "double", "double", "double", "categorical", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

      % Specify file level properties
      opts.ExtraColumnsRule = "ignore";
      opts.EmptyLineRule = "read";

      % Specify variable properties
      opts = setvaropts(opts, ["model", "cell", "radiation_type", "particle"], "EmptyFieldRule", "auto");
      opts = setvaropts(opts, "calculusType", "TrimNonNumeric", true);
      opts = setvaropts(opts, "calculusType", "ThousandsSeparator", ",");

      % Import the data
      KMKLUTLQparametersMKM = readtable(filename, opts);

      end
      
   end
end