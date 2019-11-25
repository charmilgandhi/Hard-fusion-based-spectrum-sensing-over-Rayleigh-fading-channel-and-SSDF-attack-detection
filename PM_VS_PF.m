%Comparision for Analytical and Simulational Results
%for [Pm vs Pf]

clc
close all
clear all
L = 1000;
N = 10000;
snr = -8;
snrLinear = 10.^(snr./10);
Pf = 0.01:0.01:0.5;
threshold  = (qfuncinv(Pf)./sqrt(L))+ 1;
%threshold = (qfuncinv(Pf));

for ii = 1:length(Pf)
    ii
    count = 0;
    for kk=1:10000
         n = randn(1,L);
         s = sqrt(snrLinear).*randn(1,L);
         y = s + n;
         energy = abs(y).^2; 
         energy1 =(1/L).*sum(energy);
         
         if(energy1 >= threshold(ii)) 
             count = count+1;
         end
    end
    pm_sim(ii) = count/N;
    pm_sim(ii) = 1-pm_sim(ii);
end

plot(Pf, pm_sim,'-kp','LineWidth',1)
hold on
xlabel("Prob. of false alarm (Pf)");
ylabel("Prob. of mis detection (Pm)");
title("Pm Vs. Pf");
pm_a = 1-qfunc(((threshold - (snrLinear + 1)).*sqrt(L))./(sqrt(2).*(snrLinear + 1)));
plot(Pf, pm_a,'-x','LineWidth',1.5)
legend("Simulation", "Analytical");
hold on
