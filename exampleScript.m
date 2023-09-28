matRad_rc;
matRad_cfg = MatRad_Config.instance();
matRad_cfg.propOpt.defaultMaxIter = 10000;

load('C:\r408i_data\r408i_data\CTDatasetMotion\102_HM10395_333.mat');

ct.numOfCtScen = 1;
ct.cubeHU = ct.cubeHU(1:1);
for k=1:size(cst,1)
    cst{k,4} = cst{k,4}(1:1);
end

% Lungs
cst{2,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(500,0));
cst{3,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(500,0));

% Heart

cst{4,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(500,0));
cst{13,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,0));
%% meta information for treatment plan (1) 
pln.numOfFractions  = 30;
pln.radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln.machine         = 'Generic'; %'HITfixedBL';

% beam geometry settings
pln.propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.gantryAngles    = [90]; % [?] ;
pln.propStf.couchAngles     = zeros(numel(pln.propStf.gantryAngles),1); % [?] ; 
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln.propDoseCalc.calcLET = 0;

pln.propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln.propOpt.spatioTemp      = 0;
pln.propOpt.STscenarios     = 1;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

quantityOpt  = 'physicalDose';     % options: physicalDose, effect, RBExD
%=======================================> Model check error in bioModel
modelName    = 'none';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'wcScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation

pln.multScen = matRad_multScen(ct,scenGenType);

%% stf
stf = matRad_generateStf(ct,cst,pln);

%% Dose calc
pln.propDoseCalc.clearVoxelsForRobustness = 'objectivesOnly'; % none, targetOnly, oarOnly, objectivesOnly, [scenario indexes];

pln.propDoseCalc.precalcProbabilisticQuantitites = true;
dij  = matRad_calcParticleDose(ct,stf, pln,cst,0);


%% Fluence optimization
for voiIdx=1:size(cst,1)
    if isequal(cst{voiIdx,3}, 'TARGET')
        %if ~isempty(cst{voiIdx,6})
            for objIdx=1:size(cst{voiIdx,6},2)
                
                if isa(cst{voiIdx,6}{objIdx}, 'OmegaObjectives.matRad_OmegaObjective')
                    cst{voiIdx,6}(objIdx) = [];
                else 
                    cst{voiIdx,6}{objIdx}.robustness = 'none';
                end
            end
        %end
    else
        for objIdx=1:size(cst{voiIdx,6},2)
 
                if isa(cst{voiIdx,6}{objIdx}, 'OmegaObjectives.matRad_OmegaObjective')
                    cst{voiIdx,6}(objIdx) = [];
                else 
                    cst{voiIdx,6}{objIdx}.robustness = 'none';
                end
        end
        
    end
end

tic;
resultGUI_nominal = matRad_fluenceOptimization(dij,cst,pln);
nominal_time = toc;
%% Robust planning on the scenarios
for voiIdx=1:size(cst,1)
    if isequal(cst{voiIdx,3}, 'TARGET')
        %if ~isempty(cst{voiIdx,6})
            for objIdx=1:size(cst{voiIdx,6},2)

                if isa(cst{voiIdx,6}{objIdx}, 'OmegaObjectives.matRad_OmegaObjective')
                    cst{voiIdx,6}(objIdx) = [];
                else 
                    cst{voiIdx,6}{objIdx}.robustness = 'STOCH';
                end
            end
        %end
    else
        for objIdx=1:size(cst{voiIdx,6},2)
        	     
                if isa(cst{voiIdx,6}{objIdx}, 'OmegaObjectives.matRad_OmegaObjective')
                    cst{voiIdx,6}(objIdx) = [];

                else
                    cst{voiIdx,6}{objIdx}.robustness = 'STOCH'; 
                end
        end
        
    end

end


tic;
resultGUI_stoch = matRad_fluenceOptimization(dij,cst,pln);
stoch_time = toc;
%% Robustness using Omega
for voiIdx=1:size(cst,1)
    if isequal(cst{voiIdx,3}, 'TARGET')
        %if ~isempty(cst{voiIdx,6})

        for objIdx=1:size(cst{voiIdx,6},2)
                cst{voiIdx,6}{objIdx}.robustness = 'PROB';
            end
        %end
    else
        for objIdx=1:size(cst{voiIdx,6},2)

                cst{voiIdx,6}{objIdx}.robustness = 'PROB'; 

        end
        
    end

end
for voiIdx = [2,3,4,8,13]
    cst{voiIdx,6}{2} = OmegaObjectives.matRad_TotalVariance();
    cst{voiIdx,6}{2}.penalty = 500*cst{voiIdx,6}{1}.penalty;
end
% cst{8,6}{2} = OmegaObjectives.matRad_TotalVariance();
% cst{8,6}{2}.penalty = cst{voiIdx,6}{1}.penalty;

tic;

resultGUI_prob = matRad_fluenceOptimization(dij,cst,pln);
prob_time = toc;
%% Plan comparison
clear nominalDoseDistOnScenarios;
clear probDoseDistOnScenarios;
clear stochDoseDistOnScenarios;

slice = 63;

colorCodes = [0 0.4470 0.7410; ...
    0.8500 0.3250 0.0980; ...
    0.9290 0.6940 0.1250; ...
    0.4940 0.1840 0.5560; ...
    0.4660 0.6740 0.1880; ...
    0.4660 0.6740 0.1880; ...
    0.3010 0.7450 0.9330; ...
    0.3010 0.7450 0.9330; ...
    0.6350 0.0780 0.1840; ...
    0.2010 0.7450 0.9330; ...
    0.1010 0.7450 0.9330; ...
    0.0010 0.7450 0.9330;
    0.0010 0.7450 0.7330];


doseDist = resultGUI_nominal.physicalDose;
includeStruct = [2,3,4,8,13];

lege = [];

figure;
%subplot(1,2,1);
%matRad_plotSliceWrapper(gca(),ct,cst,1,doseDist,3,slice);

dvh_nominal = matRad_calcDVH(cst,doseDist);

%subplot(1,2,2);
for structIdx=includeStruct
    plot(dvh_nominal(structIdx).doseGrid, dvh_nominal(structIdx).volumePoints, '-', 'Color', colorCodes(structIdx,:));
    hold on;
    lege = [lege, cst(structIdx,2)];
end
grid on;

grid minor;
legend(lege);
xlabel('Dose [Gy]', 'FontSize',37);
ylabel('Volume [%]', 'FontSize',37);
xlim([0,2]);
sgtitle('Nominal Plan','FontSize',37);

%%% Visualize multiple dvhs

nonEmptyScenarios = find(~cellfun(@isempty, dij.physicalDose));
w = resultGUI_nominal.w;
for scenIdx=2:numel(nonEmptyScenarios)
    tmpResultGUI = matRad_calcCubes(w, dij,nonEmptyScenarios(scenIdx));
    nominalDoseDistOnScenarios{scenIdx} = tmpResultGUI.physicalDose;%(dij.physicalDose{nonEmptyScenarios(scenIdx)}*resultGUI_nominal.w);
    dvh_nominalOnScenarios{scenIdx} = matRad_calcDVH(cst,nominalDoseDistOnScenarios{scenIdx});
end
%%% Vis
hold on;
for scenIdx=2:numel(nonEmptyScenarios)
    for structIdx=includeStruct
        plot(dvh_nominalOnScenarios{scenIdx}(structIdx).doseGrid, dvh_nominalOnScenarios{scenIdx}(structIdx).volumePoints, '--', 'Color', colorCodes(structIdx, :));
    end
end
legend(lege, 'FontSize',37);

%%% Same with robust stochastic weights
doseDist = resultGUI_stoch.physicalDose;
figure;
%subplot(1,2,1);
%matRad_plotSliceWrapper(gca(),ct,cst,1,doseDist,3,slice);

dvh_stoch = matRad_calcDVH(cst,doseDist);

%subplot(1,2,2);
for structIdx=includeStruct
    plot(dvh_stoch(structIdx).doseGrid, dvh_stoch(structIdx).volumePoints, '-', 'Color', colorCodes(structIdx,:));
    hold on;
end
grid on;

grid minor;
legend(lege);
xlabel('Dose [Gy]');
ylabel('Volume [%]');
xlim([0,2]);
sgtitle('STOCH Plan');

%%% Compute ion scenarios
nonEmptyScenarios = find(~cellfun(@isempty, dij.physicalDose));
w = resultGUI_stoch.w;
for scenIdx=2:numel(nonEmptyScenarios)
    tmpResultGUI = matRad_calcCubes(w, dij,nonEmptyScenarios(scenIdx));
    stochDoseDistOnScenarios{scenIdx} = tmpResultGUI.physicalDose;%(dij.physicalDose{nonEmptyScenarios(scenIdx)}*resultGUI_nominal.w);
    dvh_stochOnScenarios{scenIdx} = matRad_calcDVH(cst,stochDoseDistOnScenarios{scenIdx});
end
%%% Vis
hold on;

for scenIdx=2:numel(nonEmptyScenarios)
    for structIdx=includeStruct

        plot(dvh_stochOnScenarios{scenIdx}(structIdx).doseGrid, dvh_stochOnScenarios{scenIdx}(structIdx).volumePoints, '--', 'Color', colorCodes(structIdx, :));
    end
end
legend(lege);


%%% Same with OMega matrix opt
doseDist = resultGUI_prob.physicalDose;
figure;
%subplot(1,2,1);

%matRad_plotSliceWrapper(gca(),ct,cst,1,doseDist,3,slice);


dvh_prob = matRad_calcDVH(cst,doseDist);

%subplot(1,2,2);
for structIdx=includeStruct
    plot(dvh_prob(structIdx).doseGrid, dvh_prob(structIdx).volumePoints, '-', 'Color', colorCodes(structIdx,:));
    hold on;
end
grid on;

grid minor;
legend(lege);
xlabel('Dose [Gy]','FontSize',37);
ylabel('Volume [%]', 'FontSize',37);
xlim([0,2]);
sgtitle('PROB Plan', 'FontSize',37);

%%% Compute ion scenarios
nonEmptyScenarios = find(~cellfun(@isempty, dij.physicalDose));
w = resultGUI_prob.w;
for scenIdx=2:numel(nonEmptyScenarios)
    tmpResultGUI = matRad_calcCubes(w, dij,nonEmptyScenarios(scenIdx));
    probDoseDistOnScenarios{scenIdx} = tmpResultGUI.physicalDose;%(dij.physicalDose{nonEmptyScenarios(scenIdx)}*resultGUI_nominal.w);
    prob_stochOnScenarios{scenIdx} = matRad_calcDVH(cst,probDoseDistOnScenarios{scenIdx});
end

%%% Vis

hold on;
for scenIdx=2:numel(nonEmptyScenarios)
    for structIdx=includeStruct

        plot(prob_stochOnScenarios{scenIdx}(structIdx).doseGrid, prob_stochOnScenarios{scenIdx}(structIdx).volumePoints, '--', 'Color', colorCodes(structIdx, :));
    end
end
legend(lege, 'FontSize',37);

%% Compare nominal dvhs
figure;
nCols = min(3,numel(includeStruct));
nRows = ceil(numel(includeStruct)/nCols);
structNum = 1;
for structIdx=includeStruct
    subplot(nRows,nCols,structNum);
    structNum = structNum +1;
    plot(dvh_nominal(structIdx).doseGrid, dvh_nominal(structIdx).volumePoints, '-');
    hold on;
    plot(dvh_stoch(structIdx).doseGrid, dvh_stoch(structIdx).volumePoints, '-')
    plot(dvh_prob(structIdx).doseGrid, dvh_prob(structIdx).volumePoints, '-');
    legend('Nominal plan', 'Stoch plan', 'Prob plan');
    
    grid on;
    grid minor;
    xlim([0,2]);
    xlabel('Dose [Gy]');
    ylabel('Volume [%]');
    title(cst{structIdx,2});
end

sgtitle('DVH comparison');


%% Ct GIF
f = figure('WindowState', 'normal');
colormap(gray);
for k=1:10
    img = ct.cubeHU{k}(:,:,slice);
    imagesc(img);
    matRad_plotVoiContourSlice(gca(f),cst,ct,k,1,3,slice);
    %set(gcf,'MenuBar','none');
    set(gca,'DataAspectRatioMode','auto');
    set(gca,'Position',[0 0 1 1]);
    drawnow;
    fram(k) = getframe(gcf);
    im{k} = frame2im(fram(k));
end

for k=1:10
    [img, map] = rgb2ind(im{k}, 256);
    if k == 1
        imwrite(img,map,'testMethod2.gif',"gif","LoopCount",Inf,"DelayTime",0.2);
    else
        imwrite(img,map,'testMethod2.gif',"gif","WriteMode","append","DelayTime",0.2);
    end
end

%% nominal dose GIF
f = figure('WindowState', 'normal');
colormap(gray);
count = 1;
ctIdx = 1;
for k=1:size(nominalDoseDistOnScenarios,2)
    if k ==1
        matRad_plotSliceWrapper(gca(),ct,cst,1,resultGUI_nominal.physicalDose,3,slice);
        count = count +1;
    else
        if count>9
            count = 1;
            ctIdx = ctIdx+1;
        else
            count = count+1;
        end
        matRad_plotSliceWrapper(gca(),ct,cst,ctIdx,nominalDoseDistOnScenarios{k},3,slice);
    end
    %imagesc(img);
    %matRad_plotVoiContourSlice(gca(f),cst,ct,k,1,3,slice);
    set(gcf,'MenuBar','none');
    set(gca,'DataAspectRatioMode','auto');
    set(gca,'Position',[0 0 1 1]);
    drawnow;
    fram(k) = getframe(gcf);
    im{k} = frame2im(fram(k));
end

for k=1:size(nominalDoseDistOnScenarios,2)
    [img, map] = rgb2ind(im{k}, 256);
    if k == 1
        imwrite(img,map,'gifs/fsdfs.gif',"gif","LoopCount",Inf,"DelayTime",0.1);
    else
        imwrite(img,map,'gifs/fsdfs.gif',"gif","WriteMode","append","DelayTime",0.1);
    end
end

%% robust dose GIF
f = figure('WindowState', 'normal');
colormap(gray);
count = 1;
ctIdx = 1;
for k=1:size(nominalDoseDistOnScenarios,2)
    % if k ==1
    %     matRad_plotSliceWrapper(gca(),ct,cst,1,resultGUI_prob.physicalDose,3,slice);
    % else
    % 
    %     matRad_plotSliceWrapper(gca(),ct,cst,1,probDoseDistOnScenarios{k},3,slice);
    % end
     if k ==1
        matRad_plotSliceWrapper(gca(),ct,cst,1,resultGUI_nominal.physicalDose,3,slice);
        count = count +1;
    else
        if count>9
            count = 1;
            ctIdx = ctIdx+1;
        else
            count = count+1;
        end
        matRad_plotSliceWrapper(gca(),ct,cst,ctIdx,nominalDoseDistOnScenarios{k},3,slice);
    end
    set(gcf,'MenuBar','none');
    set(gca,'DataAspectRatioMode','auto');
    set(gca,'Position',[0 0 1 1]);
    drawnow;
    fram(k) = getframe(gcf);
    im{k} = frame2im(fram(k));
end

for k=1:10
    [img, map] = rgb2ind(im{k}, 256);
    if k == 1
        imwrite(img,map,'gifs/asdf.gif',"gif","LoopCount",Inf,"DelayTime",0.1);
    else
        imwrite(img,map,'gifs/asdf.gif',"gif","WriteMode","append","DelayTime",0.1);
    end
end

%% Variance maps
clear nominalScenarioDoses;
clear probScenarioDoses;
nominalScenarioDoses(:,:,:,1) = resultGUI_nominal.physicalDose;
for k =2:size(nominalDoseDistOnScenarios,2)

    nominalScenarioDoses(:,:,:,k) = nominalDoseDistOnScenarios{k};
end

nominalStd = std(nominalScenarioDoses,0,4);
figure;
matRad_plotSliceWrapper(gca(),ct,cst,1,nominalStd,3,slice);


% robust variance
probScenarioDoses(:,:,:,1) = resultGUI_prob.physicalDose;
for k =2:size(probDoseDistOnScenarios,2)
    probScenarioDoses(:,:,:,k) = probDoseDistOnScenarios{k};
end

probStd = std(probScenarioDoses,0,4);
figure;
matRad_plotSliceWrapper(gca(),ct,cst,1,probStd,3,slice);

figure;
subplot(1,2,1);
matRad_plotSliceWrapper(gca(),ct,cst,1,nominalStd,3,slice);
subplot(1,2,2);
matRad_plotSliceWrapper(gca(),ct,cst,1,probStd,3,slice);
%% GIF
%Method 1
figure;
for k=1:10
    imagesc(ct.cubeHU{k}(:,:,slice));
    exportgraphics(gcf,'testMethod1.gif', 'Append', true);
end

%% Method 2
f=figure;
map = colormap('gray');
for k=1:10
    img = (ct.cubeHU{k}(:,:,slice) - min(ct.cubeHU{k}(:,:,slice), [], 'all'));
    img = 256*img./max(img,[], 'all');

    if k == 1
        imwrite(img,map,'testMethod2.gif',"gif","LoopCount",Inf,"DelayTime",0.2, 'Screensize', size(img));
    else
        imwrite(img,map,'testMethod2.gif',"gif","WriteMode","append","DelayTime",0.2);
    end
end
%close all;
%% Method 3


f = figure('WindowState', 'normal');
colormap(gray);
for k=1:10
    img = ct.cubeHU{k}(:,:,slice);
    imagesc(img);
    matRad_plotVoiContourSlice(gca(f),cst,ct,k,1,3,slice);
    set(gcf,'MenuBar','none');
    set(gca,'DataAspectRatioMode','auto');
    set(gca,'Position',[0 0 1 1]);
    drawnow;
    fram(k) = getframe(gcf);
    im{k} = frame2im(fram(k));
end

for k=1:10
    [img, map] = rgb2ind(im{k}, 256);
    if k == 1
        imwrite(img,map,'testMethod2.gif',"gif","LoopCount",Inf,"DelayTime",0.2);
    else
        imwrite(img,map,'testMethod2.gif',"gif","WriteMode","append","DelayTime",0.2);
    end
end


% writerObj = VideoWriter('testMetghod3.avi');
% writerObj.FrameRate = 10;
% open(writerObj);
% for k=1:length(fram)
%     writeVideo(writerObj, fram(k));
% end
% close(writerObj);

%% asdfa
for i = 1:N
    figure(1)  
    imshow(processo(:,:,1,i))
      hold on
      plot(X,Y,'o')
      plot(X0,Y0,'o')
      plot(X1,Y1,'o')
      plot(X2,Y2,'o')
      plot(X3,Y3,'o')
      hold off
      F(i) = getframe(gcf) ;
      drawnow
    end
  % create the video writer with 1 fps
  % set the seconds per image
% open the video writer
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);
