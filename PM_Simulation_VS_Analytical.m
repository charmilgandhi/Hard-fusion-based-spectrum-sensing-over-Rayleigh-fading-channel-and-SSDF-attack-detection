%Comparision of Simulational and Analytical results
% of Pm vs SNR

clc
close all
clear all
L = 50;
monte_carlo =10^5;
Pf = 0.2;
snr_db = -10:1:20;
snr = 10.^(snr_db./10); 
threshold = 2*gammaincinv(1-Pf,L/2);
n=randn(monte_carlo,L);

%Simulated probability of detection for Rayleigh channel
for tt = 1:length(snr_db)
    tt
    s = [];
    h = []; 
    mes = randi([0 1],monte_carlo,L);
    s = (2.*(mes)-1); 
    h1 = (randn(monte_carlo,1)+1j*randn(monte_carlo,1))./(sqrt(2)); 
    h = repmat(h1,1,L);
    y = sqrt(snr(tt)).*abs(h).*s + n;
    energy = (abs(y).^2); 
    for kk=1:monte_carlo 
        energy_final(kk) =sum(energy(kk,:));
    end
    temp = (energy_final >= threshold);
    i = sum(temp);
    Pd_sim(tt) = i/kk; 
    Pm_sim(tt) = 1-Pd_sim(tt);
end

plot(snr_db,Pm_sim,'-o','LineWidth',2,'color','red')
hold on
grid on
axis([-10 10 0 1]);

Pd_theory  = [];
u = L./2; % Time-Bandwidth product

for ii = 1:length(snr_db)
    temp1 = (L/2)*snr(ii);
    A = temp1./(1 + temp1);
    n = 0:1:u-2;
    t1 = sum((1./factorial(n)).*(threshold./2).^(n));
    t2 = sum((1./factorial(n)).*(((threshold./2).*(A)).^(n)));
    Pd_theory(ii) = exp(-threshold./2).*t1 + (1./A).^(u-1).*(exp(-threshold./(2.*(1+temp1))) - exp(-threshold./2).*t2);
    pm_theory(ii) = 1-Pd_theory(ii);
end

plot(snr_db,pm_theory,'-x','LineWidth',1,'color','Green');
legend('Pm (Simulation)','Pm (Theory)')
xlabel('SNR (dB)')
ylabel('Probability Of Miss Detection(Pm)')
title('Pm vs SNR');