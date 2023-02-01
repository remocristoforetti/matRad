track = Track_KC();
track.Z_ion = 1;

domain = Domain();
domain.rNuc = 0.32;

InteGral = IntegralDose();
InteGral.Resolution = 0.01;
InteGral.stepWeightedIntegral = 0.01;
InteGral.corFact = 1;

E = logspace(-1,log10(2*10^(2)), 70);
g = [];
for k=1:size(E,2)
    track.Ek = E(k)*track.A;

    InteGral.computeZD(track,domain);
    g(k) = InteGral.gamma;
    k
end
EA = logspace(-2,-1,70);
gA = [];
for k=1:size(EA,2)
    track.Ek = EA(k)*track.A;
    InteGral.computeZD(track,domain);
    
    gA(k) = InteGral.gamma;
    k
end

ET = [EA,E];

GT = [gA,g];

figure;

semilogx(ET,GT,'.-');
grid on;