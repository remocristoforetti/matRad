matRad_rc;
load('TG119.mat');

pln.radiationMode = 'protons';        
pln.machine       = 'generic_MCsquare';

pln.propOpt.bioOptimization = 'physicalDose';

pln.propDoseCalc.engine = 'FRED';

pln.propDoseCalc.calcLET = 0;


pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90,90];
pln.propStf.couchAngles   = [0,45];
pln.propStf.bixelWidth    = 3;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

%testStfSingleField;

stf = matRad_generateStf(ct,cst,pln);

%pln.propDoseCalc.engine = 'HongPB';
%dij = matRad_calcDoseInfluence(ct,cst,stf,pln);

%resultGUI = matRad_fluenceOptimization(dij,cst,pln);
% w = resultGUI.w;
% save('w.mat', 'w');

% % %% add random weights
%  for i=1:numel(stf)
%      for j=1:stf(i).numOfRays
%          stf(i).ray(j).weight = 3*rand(1,stf(i).numOfBixelsPerRay(j),'double');
%      end
%  end

% load('w.mat');
% 
% 
% counter = 0;
% counterBixel=1;
% for i=1:numel(stf)
%     for j=1:stf(i).numOfRays
%         stf(i).ray(j).weight = w(counter+1:counter+stf(i).numOfBixelsPerRay(counterBixel));%3*rand(1,stf(i).numOfBixelsPerRay(j),'double');
%         counter = counter+stf(i).numOfBixelsPerRay(counterBixel);
%         counterBixel = counterBixel+1;
%     end
% end
% w = [];
% for i=1:size(stf,1)
%     for j =1:stf(i).numOfRays 
%         w = [w, stf(i).ray(j).weight];
%     end
% end

pln.propDoseCalc.engine = 'HongPB';

% matRad_cfg = MatRad_Config.instance();
% matRad_cfg.propDoseCalc.defaultNumHistoriesDirect = 1000;
% matRad_cfg.propDoseCalc.defaultNumHistoriesPerBeamlet = 1000;
dij_PB = matRad_calcDoseInfluence(ct,cst,stf,pln);

resultGUI = matRad_fluenceOptimization(dij_PB,cst,pln);

% %% Dose calc MCsquare
% matRad_cfg = MatRad_Config.instance();
% 
% matRad_cfg.propDoseCalc.defaultNumHistoriesPerBeamlet = 10000*dij_PB.totalNumOfBixels;%floor((1000*dij_PB.totalNumOfBixels)/sum(resultGUI.w));
% 
% pln.propDoseCalc.engine = 'MCsquare';
% 
% resultGUI_MCsquare = matRad_calcDoseDirectMC(ct,stf,pln,cst,resultGUI.w);
%% Dose Calculation
stfFred = stf;
counter = 0;


for i=1:numel(stfFred)
    counterBixel=1;
    for j=1:stfFred(i).numOfRays
        
        stfFred(i).ray(j).weight = resultGUI.w(counter+1:counter+stfFred(i).numOfBixelsPerRay(counterBixel));%3*rand(1,stf(i).numOfBixelsPerRay(j),'double');
        % if i==1
        %      stfFred(i).ray(j).weight = 0*stfFred(i).ray(j).weight;
        % end
        counter = counter+stfFred(i).numOfBixelsPerRay(counterBixel);
        counterBixel = counterBixel+1;
    end
end
pln.propDoseCalc.engine = 'FRED';

dij = matRad_calcDoseInfluence(ct,cst,stfFred,pln);
%% read


cube = matRad_readMhd('FRED/MCrun/out', 'Dose.mhd');

cube = matRad_interp3(dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z, ...
                                             cube, ...
                                             dij.ctGrid.x,dij.ctGrid.y',dij.ctGrid.z,'linear',0);


cubePB = reshape(dij_PB.physicalDose{1}*resultGUI.w, dij_PB.doseGrid.dimensions);
cubePB = matRad_interp3(dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z, ...
                                             cubePB, ...
                                             dij.ctGrid.x,dij.ctGrid.y',dij.ctGrid.z,'linear',0);
%% Visualize
slice1 = 83;
slice2 = 83;
slice3 = 53;

figure;
movegui(gcf(),'southwest');
subplot(3,3,1);
matRad_plotSliceWrapper(gca(),ct,cst,1,cube,1,slice1);

subplot(3,3,2);
matRad_plotSliceWrapper(gca(),ct,cst,1,cube,2,slice2);%
%imagesc(squeeze(cube(:,83,:)));
%matRad_plotVoiContourSlice(gca(),cst,ct.cubeHU,1,[],2,83);

subplot(3,3,3);
matRad_plotSliceWrapper(gca(),ct,cst,1,cube,3,slice3);


subplot(3,3,4);
matRad_plotSliceWrapper(gca(),ct,cst,1,cubePB,1,slice1);

subplot(3,3,5);
matRad_plotSliceWrapper(gca(),ct,cst,1,cubePB,2,slice2);%
%imagesc(squeeze(cube(:,83,:)));
%matRad_plotVoiContourSlice(gca(),cst,ct.cubeHU,1,[],2,83);

subplot(3,3,6);
matRad_plotSliceWrapper(gca(),ct,cst,1,cubePB,3,slice3);

gamma1 = matRad_gammaIndex(cube,cubePB,[ct.resolution.x,ct.resolution.y,ct.resolution.z],[3 3]);
%cMaApGamma = colormap(gammaIndex);
subplot(3,3,7);
matRad_plotSliceWrapper(gca(),ct,cst,1,gamma1,1,slice1,[],[],[],gammaIndex,[0,2]);

subplot(3,3,8);
matRad_plotSliceWrapper(gca(),ct,cst,1,gamma1,2,slice2,[],[],[],gammaIndex,[0,2]);

subplot(3,3,9);
matRad_plotSliceWrapper(gca(),ct,cst,1,gamma1,3,slice3,[],[],[],gammaIndex,[0,2]);

%% Visualize
slice1 = 83;
slice2 = 83;
slice3 = 53;

figure;
movegui(gcf(),'southwest');
subplot(3,3,1);
matRad_plotSliceWrapper(gca(),ct,cst,1,cube,1,slice1);

subplot(3,3,2);
matRad_plotSliceWrapper(gca(),ct,cst,1,cube,2,slice2);%
%imagesc(squeeze(cube(:,83,:)));
%matRad_plotVoiContourSlice(gca(),cst,ct.cubeHU,1,[],2,83);

subplot(3,3,3);
matRad_plotSliceWrapper(gca(),ct,cst,1,cube,3,slice3);


subplot(3,3,4);
matRad_plotSliceWrapper(gca(),ct,cst,1,resultGUI_MCsquare.physicalDose,1,slice1);

subplot(3,3,5);
matRad_plotSliceWrapper(gca(),ct,cst,1,resultGUI_MCsquare.physicalDose,2,slice2);%
%imagesc(squeeze(cube(:,83,:)));
%matRad_plotVoiContourSlice(gca(),cst,ct.cubeHU,1,[],2,83);

subplot(3,3,6);
matRad_plotSliceWrapper(gca(),ct,cst,1,resultGUI_MCsquare.physicalDose,3,slice3);

gamma1 = matRad_gammaIndex(cube,resultGUI_MCsquare.physicalDose,[ct.resolution.x,ct.resolution.y,ct.resolution.z],[3 3]);
%cMaApGamma = colormap(gammaIndex);
subplot(3,3,7);
matRad_plotSliceWrapper(gca(),ct,cst,1,gamma1,1,slice1,[],[],[],gammaIndex,[0,2]);

subplot(3,3,8);
matRad_plotSliceWrapper(gca(),ct,cst,1,gamma1,2,slice2,[],[],[],gammaIndex,[0,2]);



subplot(3,3,9);
matRad_plotSliceWrapper(gca(),ct,cst,1,gamma1,3,slice3,[],[],[],gammaIndex,[0,2]);