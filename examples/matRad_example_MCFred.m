matRad_rc;
%load('TG119.mat');
load('BOXPHANTOM.mat');

pln.radiationMode = 'protons';

pln.machine       = 'generic_MCsquare';
pln.propOpt.bioOptimization = 'physicalDose';

pln.propDoseCalc.calcLET = 0;

pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [0];
pln.propStf.couchAngles   = [0];
pln.propStf.bixelWidth    = 3;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 2; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 2; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 2; % [mm]


%testStfSingleField;

%%% reduce target to get less bixels
% maskV = zeros(ct.cubeDim);
% targetIndex = find(strcmp([cst(:,3)], 'TARGET'));
% maskV(cst{targetIndex,4}{1}) = 1;
% 
% c = strel('disk',10);
% maskV = imerode(maskV,c);
% 
% cst{targetIndex,4}{1} = find(maskV);
% %%%
%stf = matRad_generateStf(ct,cst,pln);

machine = matRad_loadMachine(pln);
%stf = matRad_generateStfSinglePencilBeam(ct,cst,pln, machine.data(37).energy, 0, 0);
x = [-100:30:100];
y = [-100:30:100];


stf = matRad_generateStfSpotGridForTesting(ct,cst,pln,machine.data(37).energy, x,y);


%pln.bioParam = matRad_bioModel(pln.radiationMode, );
pln.propDoseCalc.engine = 'HongPB';
pln.propDoseCalc.calcLET = false;
dij_PB = matRad_calcDoseInfluence(ct,cst,stf,pln);

%resultGUI = matRad_fluenceOptimization(dij_PB,cst,pln);

% %% Dose calc MCsquare
matRad_cfg = MatRad_Config.instance();


matRad_cfg.propDoseCalc.defaultNumHistoriesPerBeamlet = 1000;%1000*dij_PB.totalNumOfBixels;%floor((1000*dij_PB.totalNumOfBixels)/sum(resultGUI.w));

pln.propDoseCalc.engine = 'MCsquare';


%dij_MCsquare = matRad_calcDoseInfluence(ct,cst,stf,pln);

resultGUI_MCsquare = matRad_calcDoseDirect(ct,stf,pln,cst,ones(stf.totalNumOfBixels,1));
%% Dose Calculation
stfFred = stf;
counter = 0;
calcDoseDirect = 1;

% if calcDoseDirect == 1
% 
%     for i=1:numel(stfFred)
% 
%         counterBixel=1;
%         for j=1:stfFred(i).numOfRays
% 
%             stfFred(i).ray(j).weight = resultGUI.w(counter+1:counter+stfFred(i).numOfBixelsPerRay(counterBixel));%3*rand(1,stf(i).numOfBixelsPerRay(j),'double');
%             % if i==1
%             %      stfFred(i).ray(j).weight = 0*stfFred(i).ray(j).weight;
%             % end
%             counter = counter+stfFred(i).numOfBixelsPerRay(counterBixel);
% 
%             counterBixel = counterBixel+1;
% 
%         end
%     end
% end
pln.propDoseCalc.engine = 'FRED';
pln.propDoseCalc.calcLET = true;

resultGUI.w = ones(stf.totalNumOfBixels,1);


pln.propOpt.bioOptimization = 'MCN_RBExD';
resultGUI = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);
%dij = matRad_calcDoseInfluence(ct,cst,stfFred,pln);

%resultGUI.w = 1;
%% read

%cube = matRad_readMhd('FRED/MCrun/out', 'Dose.mhd');

%cube = matRad_interp3(dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z, ...
%                                             cube, ...
%                                             dij.ctGrid.x,dij.ctGrid.y',dij.ctGrid.z,'linear',0);

try
    cube = reshape(full(dij.physicalDose{1}*resultGUI.w),dij.doseGrid.dimensions);
    %cube = matRad_interp3(dij.doseGrid.x, dij.doseGrid.y', dij.doseGrid.z, cube, dij.ctGrid.x, dij.ctGrid.y', dij.ctGrid.z, 'linear',0);
    
    cubePB = reshape(full(dij_PB.physicalDose{1}*resultGUI.w), dij_PB.doseGrid.dimensions);
    
    %cubePB = matRad_interp3(dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z, ...
    %                                             cubePB, ...
    %                                             dij.ctGrid.x,dij.ctGrid.y',dij.ctGrid.z,'linear',0);
    
    cube_MCsquare = reshape(full(dij_MCsquare.physicalDose{1}*resultGUI.w), dij_MCsquare.doseGrid.dimensions);
    
    %cube_MCsquare = matRad_interp3(dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z, ...
    %                                             cube_MCsquare, ...
    %                                             dij.ctGrid.x,dij.ctGrid.y',dij.ctGrid.z,'linear',0);
    currCT = ct;
    currCT.cube{1} = matRad_interp3(dij.ctGrid.x,dij.ctGrid.y',dij.ctGrid.z, ...
                                                 ct.cube{1}, ...
                                                 dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'linear',0);
    currCT.cubeHU{1} = matRad_interp3(dij.ctGrid.x,dij.ctGrid.y',dij.ctGrid.z, ...
                                                 ct.cubeHU{1}, ...
                                                 dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'linear',0);
    
    currCT.cubeDim = size(cube);
    currCT.resolution = dij.doseGrid.resolution;
    
    cstOnGrid = matRad_resizeCstToGrid(cst,dij.ctGrid.x,dij.ctGrid.y, dij.ctGrid.z, dij.doseGrid.x, dij.doseGrid.y, dij.doseGrid.z);
catch
    cube = resultGUI.physicalDose;
    cubeLET = resultGUI.LET;
    %cubeLET_PB = dij_PB.mLETDose{1}*resultGUI.w;
    cubePB = reshape(full(dij_PB.physicalDose{1}*resultGUI.w), dij_PB.doseGrid.dimensions);
    
    cubePB = matRad_interp3(dij_PB.doseGrid.x,dij_PB.doseGrid.y',dij_PB.doseGrid.z, ...
                                                 cubePB, ...
                                                 dij_PB.ctGrid.x,dij_PB.ctGrid.y',dij_PB.ctGrid.z,'linear',0);

    % cube_MCsquare = reshape(full(dij_PB.physicalDose{1}*resultGUI.w), dij_PB.doseGrid.dimensions);
    % 
    % cube_MCsquare = matRad_interp3(dij_PB.doseGrid.x,dij_PB.doseGrid.y',dij_PB.doseGrid.z, ...
    %                                              cube_MCsquare, ...
    %                                              dij_PB.ctGrid.x,dij_PB.ctGrid.y',dij_PB.ctGrid.z,'linear',0);
    cube_MCsquare = resultGUI_MCsquare.physicalDose;

    currCT = ct;
    cstOnGrid = cst;

    cubeMCN = resultGUI.BioDose;

end
%% Get profiles
sliceShallow = 40;%93;%121;%40;%70;%121;
sliceISO = 80;%121;%250;%60;%121;%250;
profileIdx = [1:size(cube,3)];%floor(size(cubePB,2)/2);

cube1 = cube;
cube2 = cubePB;
cube3 = cubeMCN;

leg = {'FRED', 'PB', 'MCsquare'};

figure;
tiledlayout(5,2)
nexttile;

matRad_plotSliceWrapper(gca(),currCT,cstOnGrid,1,cube1,1,sliceShallow);
title(['D = ', num2str(sliceShallow)]);
ylabel('FRED', 'FontSize', 20, 'HorizontalAlignment','center');
yline(profileIdx([1,end]));

nexttile;
matRad_plotSliceWrapper(gca(),currCT,cstOnGrid,1,cube1,1,sliceISO);
title(['D = ', num2str(sliceISO)]);

nexttile;
matRad_plotSliceWrapper(gca(),currCT,cstOnGrid,1,cube2,1,sliceShallow);

title(['D = ', num2str(sliceShallow)]);
ylabel('PB', 'FontSize', 20, 'HorizontalAlignment','center');
yline(profileIdx([1,end]));

nexttile;
matRad_plotSliceWrapper(gca(),currCT,cstOnGrid,1,cube2,1,sliceISO);
title(['D = ', num2str(sliceISO)]);



%x = dij_PB.doseGrid.x - (dij_PB.doseGrid.x(end) - dij_PB.doseGrid.x(1))/2 + dij_PB.doseGrid.resolution.x;%[1:size(cube,1)]*ct.resolution.x - ct.resolution.x/2;
x = dij_PB.ctGrid.x - (dij_PB.ctGrid.x(end) - dij_PB.ctGrid.x(1))/2 + dij_PB.ctGrid.resolution.x;%[1:size(cube,1)]*ct.resolution.x - ct.resolution.x/2;

nexttile([2,1]);
profileSHALLOW = sum(squeeze(cube1(sliceShallow,profileIdx,:)),1);
profileSHALLOW2 = sum(squeeze(cube2(sliceShallow,profileIdx,:)),1);
semilogy(x, profileSHALLOW, '.-');
hold on;
semilogy(x,profileSHALLOW2, '.-');
grid on;

if exist('cube3', 'var')
    profileSHALLOW3 = sum(squeeze(cube3(sliceShallow,profileIdx,:)),1);
    semilogy(x, profileSHALLOW3, '.-');
end

legend(leg);
xlim([x(find(profileSHALLOW2>0,1,'first')-10), x(find(profileSHALLOW2>0,1,'last')+10)]);



nexttile([2,1]);
profileISO = sum(squeeze(cube(sliceISO,profileIdx,:)),1);

profileISOPB = sum(squeeze(cubePB(sliceISO,profileIdx,:)),1);

semilogy(x, profileISO, '.-');
hold on;
semilogy(x, profileISOPB, '.-');
if exist('cube3', 'var')
    profileISO3 = sum(squeeze(cube3(sliceShallow,profileIdx,:)),1);
    semilogy(x, profileISO3, '.-');
end

grid on;
legend(leg);
xlim([x(find(profileISOPB>0,1,'first')-10), x(find(profileISOPB>0,1,'last')+10)]);

nexttile;
plot(x, (profileSHALLOW - profileSHALLOW2)./profileSHALLOW2, '.-');
hold on;
plot(x, (profileSHALLOW - profileSHALLOW3)./profileSHALLOW2, '.-');
grid on;
grid minor;
legend('FRED - PB', 'FRED - MC2');
xlim([x(find(profileSHALLOW2>0,1,'first')-10), x(find(profileSHALLOW2>0,1,'last')+10)]);


nexttile;
plot(x, (profileISO - profileISOPB)./profileISOPB, '.-');
hold on;
plot(x, (profileISO - profileISO3)./profileISOPB, '.-');
grid on;
grid minor;
legend('FRED - PB', 'FRED - MC2');
xlim([x(find(profileISOPB>0,1,'first')-10), x(find(profileISOPB>0,1,'last')+10)]);

%% BP

BPprofileFred = squeeze(sum(cube1,[2,3]));
BPprofilePB = squeeze(sum(cube2,[2,3]));
BPprofileMCsquare = squeeze(sum(cube3, [2,3]));

figure;
plot(dij_PB.ctGrid.z,BPprofileFred, '.-');
hold on;
plot(dij_PB.ctGrid.z,BPprofilePB, '.-');
grid on;

plot(dij_PB.ctGrid.z,BPprofileMCsquare, '.-');

%xline(sliceShallow*currCT.resolution.y);

%xline(sliceISO*currCT.resolution.y);

legend('FRED', 'PB', 'MCsquare');
xlim([100,300]);


%% LET profiles
LETprofileFred = squeeze(sum(cubeLET,[2,3]));
figure;
plot(dij_PB.ctGrid.z,LETprofileFred, '.-');
grid on;
xlim([100,300]);



%% Visualize


slice1 = 45;
slice2 = 80;
slice3 = 80;

figure;

movegui(gcf(),'southwest');
subplot(3,3,1);
matRad_plotSliceWrapper(gca(),ct,cst,1,cube,1,slice1);
ylabel('FRED', 'FontSize', 20, 'HorizontalAlignment','center');

subplot(3,3,2);
matRad_plotSliceWrapper(gca(),ct,cst,1,cube,2,slice2);%
%imagesc(squeeze(cube(:,83,:)));
%matRad_plotVoiContourSlice(gca(),cst,ct.cubeHU,1,[],2,83);

subplot(3,3,3);
matRad_plotSliceWrapper(gca(),ct,cst,1,cube,3,slice3);


subplot(3,3,4);
matRad_plotSliceWrapper(gca(),ct,cst,1,cubePB,1,slice1);
ylabel('PB', 'FontSize', 20, 'HorizontalAlignment','center');


subplot(3,3,5);
matRad_plotSliceWrapper(gca(),ct,cst,1,cubePB,2,slice2);
%imagesc(squeeze(cube(:,83,:)));
%matRad_plotVoiContourSlice(gca(),cst,ct.cubeHU,1,[],2,83);

subplot(3,3,6);
matRad_plotSliceWrapper(gca(),ct,cst,1,cubePB,3,slice3);

gamma1 = matRad_gammaIndex(cube,cubePB,[ct.resolution.x,ct.resolution.y,ct.resolution.z],[3 3]);
%cMaApGamma = colormap(gammaIndex);
subplot(3,3,7);
matRad_plotSliceWrapper(gca(),ct,cst,1,gamma1,1,slice1,[],[],[],gammaIndex,[0,2]);
ylabel('Gamma', 'FontSize', 20, 'HorizontalAlignment','center');

subplot(3,3,8);
matRad_plotSliceWrapper(gca(),ct,cst,1,gamma1,2,slice2,[],[],[],gammaIndex,[0,2]);

subplot(3,3,9);
matRad_plotSliceWrapper(gca(),ct,cst,1,gamma1,3,slice3,[],[],[],gammaIndex,[0,2]);

sgtitle('FRED vs PB', 'FontSize', 17);

%% Visualize

slice1 = 83;
slice2 = 83;
slice3 = 53;

figure;
%movegui(gcf(),'southwest');
subplot(3,3,1);

matRad_plotSliceWrapper(gca(),ct,cst,1,cube,1,slice1);

ylabel('FRED', 'FontSize', 20, 'HorizontalAlignment','center');

subplot(3,3,2);
matRad_plotSliceWrapper(gca(),ct,cst,1,cube,2,slice2);%
%imagesc(squeeze(cube(:,83,:)));
%matRad_plotVoiContourSlice(gca(),cst,ct.cubeHU,1,[],2,83);

subplot(3,3,3);
matRad_plotSliceWrapper(gca(),ct,cst,1,cube,3,slice3);


subplot(3,3,4);
matRad_plotSliceWrapper(gca(),ct,cst,1,resultGUI_MCsquare.physicalDose,1,slice1);
ylabel('MCsquare', 'FontSize', 20, 'HorizontalAlignment','center');

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
ylabel('Gamma', 'FontSize', 20, 'HorizontalAlignment','center');

subplot(3,3,8);
matRad_plotSliceWrapper(gca(),ct,cst,1,gamma1,2,slice2,[],[],[],gammaIndex,[0,2]);



subplot(3,3,9);
matRad_plotSliceWrapper(gca(),ct,cst,1,gamma1,3,slice3,[],[],[],gammaIndex,[0,2]);


%% Dij

fNameDij = 'C:\r408i_data\r408i_data\matRad_dev\FRED\MCrun\out\Dij.bin';

fID = fopen(fNameDij, 'r');
dijTest = fread(fID);
dijTest = reshape(dijTest,[],31669);
fclose(fID);

