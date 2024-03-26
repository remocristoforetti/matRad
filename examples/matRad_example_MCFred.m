matRad_rc;
%load('TG119.mat');
load('BOXPHANTOM.mat');

pln.radiationMode = 'protons';

pln.machine       = 'newGeneric_4Aug';%'generic_MCsquare';%'newGeneric_4Aug';

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

pln.bioParam = matRad_bioModel(pln.radiationMode, 'physicalDose','none');
pln.multScen = matRad_multScen(ct, 'nomScen');

machine = matRad_loadMachine(pln);
%stf = matRad_generateStfSinglePencilBeam(ct,cst,pln, machine.data(37).energy, 0, 0);
x = [0,20,40,70]; %[-100:30:100];
y = [0]; %[-100:30:100];

stf = matRad_generateStfSpotGridForTesting(ct,cst,pln,machine.data(37).energy, x,y);


%pln.bioParam = matRad_bioModel(pln.radiationMode, );
pln.propDoseCalc.engine = 'HongPB';
pln.propDoseCalc.calcLET = false;
dij_PB = matRad_calcDoseInfluence(ct,cst,stf,pln);

%resultGUI = matRad_fluenceOptimization(dij_PB,cst,pln);

%% Dose calc MCsquare
matRad_cfg = MatRad_Config.instance();
matRad_cfg.propDoseCalc.defaultNumHistoriesDirect = 100000*sum([stf.totalNumOfBixels]);

% load([pln.radiationMode, '_', pln.machine]);
% machine.meta.created_by = machine.meta.createdBy;
% machine.meta.created_on = machine.meta.createdOn;
% save('basedata/protons_newGeneric_4Aug.mat', 'machine');

pln.propDoseCalc.engine = 'MCsquare';


dij_MCsquare = matRad_calcDoseInfluence(ct,cst,stf,pln);

%resultGUI_MCsquare = matRad_calcDoseDirect(ct,stf,pln,cst,ones(stf.totalNumOfBixels,1));
%% Dose Calculation
stfFred = stf;
counter = 0;
calcDoseDirect = 1;

pln.propDoseCalc.engine = 'FRED';
pln.propDoseCalc.calcLET = false;
matRad_cfg.propDoseCalc.defaultNumHistoriesPerBeamlet = 1000000;%*sum([stf.totalNumOfBixels]);

resultGUI.w = ones(stf.totalNumOfBixels,1);
resultGUI_FRED = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);
%dij_FRED = matRad_calcDoseInfluence(ct,cst,stf,pln);
%resultGUI = matRad_calcCubes(resultGUI.w, dij_FRED,1);
%resultGUI = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);


%% read
usedEngines = {'PB', 'FRED'};


for k=1:numel(usedEngines)

    if exist(['dij_',usedEngines{k}], 'var')
        currDij = eval(['dij_',usedEngines{k}]);
        eval(['cube', usedEngines{k}, '= reshape(full(currDij.physicalDose{1}*resultGUI.w), currDij.doseGrid.dimensions);']);
    end
end

currCT = ct;
currCT.cube{1} = matRad_interp3(dij_PB.ctGrid.x,dij_PB.ctGrid.y',dij_PB.ctGrid.z, ...
                                             ct.cube{1}, ...
                                             dij_PB.doseGrid.x,dij_PB.doseGrid.y',dij_PB.doseGrid.z,'linear',0);
currCT.cubeHU{1} = matRad_interp3(dij_PB.ctGrid.x,dij_PB.ctGrid.y',dij_PB.ctGrid.z, ...
                                             ct.cubeHU{1}, ...
                                             dij_PB.doseGrid.x,dij_PB.doseGrid.y',dij_PB.doseGrid.z,'linear',0);

currCT.cubeDim = dij_PB.doseGrid.dimensions;
currCT.resolution = dij_PB.doseGrid.resolution;


cstOnGrid = matRad_resizeCstToGrid(cst,dij_PB.ctGrid.x,dij_PB.ctGrid.y, dij_PB.ctGrid.z, dij_PB.doseGrid.x, dij_PB.doseGrid.y, dij_PB.doseGrid.z);

%% Get profiles
sliceDepth = 71;%121;%40;%121;%40;%121;%93;%121;%40;%70;%121;
sliceLateral = 121;%240;%80;%242;%80;%242;%121;%250;%60;%121;%250;
profileIdx = [1:currCT.cubeDim(3)];%floor(size(cubePB,2)/2);


figure;
%tiledlayout(numel(usedEngines)+1,2);


for k=1:numel(usedEngines)
    curCube = eval(['cube', usedEngines{k}]);
    %nexttile(k);
    subplot(numel(usedEngines)+1,2,2*k-1);
    matRad_plotSliceWrapper(gca(),currCT,cstOnGrid,1,curCube,1,sliceDepth);
    ylabel(usedEngines{k}, 'FontSize', 17);
    %nexttile([k,2]);
    subplot(numel(usedEngines)+1,2,2*k);
    matRad_plotSliceWrapper(gca(),currCT,cstOnGrid,1,curCube,3,sliceLateral);

    currProfile_1 = sum(squeeze(curCube(sliceDepth,profileIdx,:)),2);
    x_1 = [1:currCT.cubeDim(1)]*currCT.resolution.y - (currCT.cubeDim(1)*currCT.resolution.y/2);% + currCT.resolution.y;

    %nexttile([numel(usedEngines)+1,1]);
    subplot(numel(usedEngines)+1,2,2*numel(usedEngines)+1);
    hold on;
    semilogy(x_1, currProfile_1, '.-');
    grid on;
    hold off;
    xlim([-50,50])
    xlabel('mm', 'FontSize',17);
    legend(usedEngines(1:k));


    currProfile_2 = sum(squeeze(curCube(:,:,:)),[2,3]);
    x_2 = [1:currCT.cubeDim(2)]*currCT.resolution.x - currCT.resolution.x/2;
    
    %nexttile([numel(usedEngines)+1,2]);
    subplot(numel(usedEngines)+1,2,2*numel(usedEngines)+2);
    hold on;
    semilogy(x_2, currProfile_2, '.-');
    grid on;
    hold off;
    xlabel('mm', 'FontSize',17);
    legend(usedEngines(1:k));

end

%% other
figure;

tiledlayout(2,1);
nexttile;
for k=1:numel(usedEngines)
    curCube = eval(['cube', usedEngines{k}]);
    currProfile_1 = squeeze(sum(curCube(sliceDepth,profileIdx,:),2));
    x_1 = [1:currCT.cubeDim(1)]*currCT.resolution.y - (currCT.cubeDim(1)*currCT.resolution.y/2);% + currCT.resolution.y;
    hold on;
    semilogy(x_1, currProfile_1, '.-');
    grid on;
    hold off;
    xlim([-50,50])
    xlabel('mm', 'FontSize',17);
    legend(usedEngines(1:k));
end

nexttile;
for k=1:numel(usedEngines)
    
    curCube = eval(['cube', usedEngines{k}]);
    currProfile_1 = squeeze(sum(curCube(sliceDepth,profileIdx,:),3));
    
    x_1 = [1:currCT.cubeDim(1)]*currCT.resolution.y - (currCT.cubeDim(1)*currCT.resolution.y/2);% + currCT.resolution.y;
    hold on;
    semilogy(x_1, currProfile_1, '.-');
    grid on;
    hold off;
    xlim([-50,50])
    xlabel('mm', 'FontSize',17);
    legend(usedEngines(1:k));
end
%% Comp dose 


for k=1:numel(usedEngines)

    if exist(['dij_',usedEngines{k}], 'var')
        currDij = eval(['dij_',usedEngines{k}]);
        eval(['resampled', usedEngines{k}, '= matRad_calcCubes(ones(4,1),currDij,1);']);
    end
end


matRad_compareDose(resampledPB.physicalDose,resultGUI_FRED.physicalDose,ct,cst);
%% asd


cube1 = cube;
cube2 = cubePB;
cube3 = cubeMCN;

leg = {'FRED', 'PB', 'MCsquare'};

figure;
tiledlayout(5,2);
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

