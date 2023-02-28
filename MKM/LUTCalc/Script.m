track = Track_KCGPU();
Z = [1 2 3 4 5 6];
domain = Domain();
domain.rNuc = 0.35;
domain.Beta_Tissue = 0.05;

domain.RN = 3.9;

InteGral = IntegralDose();
InteGral.Resolution = 0.01;
InteGral.stepWeightedIntegral = 0.01;
InteGral.corFact = 1;

E = logspace(-1,log10(100),70);
gputime = GenerateLUT(domain, InteGral, track, 1, E, 'D:\matRad_gitHubRemo_MKM_TOPAS\MKM\LUTCalc_GPU');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

track = Track_KC();
time = GenerateLUT(domain, InteGral, track, 1, E, 'D:\matRad_gitHubRemo_MKM_TOPAS\MKM\LUTCalc');
gputime
time
%% Plot Luts
filenames = dir('X:\matRad_dev_varRBE_robOpt_TopasMKM\MKM\LUTs\*.mat');
figure;
leg = [];
for k=1:size(filenames,1)
    load(strcat(filenames(1).folder, filesep, filenames(k).name));
    leg = [leg, {strcat('Z = ', num2str(k))}];
    Lut(k).Z = k;
    Lut(k).E = E;

    Lut(k).zD = zD;
    semilogx(E,zD, '.-');
    hold on;
end

legend(leg);
grid on;

%% Write TOPAS Luts

ScorerName = {'Patient/ZMIX/'};
filename = 'IonLuts.txt';
Ions = {'protons', 'He', 'Li', 'Be', 'B', 'C'};

fid = fopen(filename, 'w');
for k=1:size(Ions,2)
    fprintf(fid, strcat('sv:Sc/', ScorerName{1}, Ions{k}, '/Energy = %i'), size(Lut(k).E,2));
    fprintf(fid, '% 3.3f ', Lut(k).E);
    fprintf(fid, '\n');
    fprintf(fid, strcat('sv:Sc/', ScorerName{1}, Ions{k}, '/zD = %i'), size(Lut(k).zD,2));
    fprintf(fid, ' %3.3f ', Lut(k).zD);
    fprintf(fid, '\n \n');
end

fclose(fid);