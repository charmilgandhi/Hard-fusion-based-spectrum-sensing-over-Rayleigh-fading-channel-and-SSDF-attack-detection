%Prob of Global Miss-Detection vs Number of SVUs

clc
close all
clear all
L = 50; 
iter =10^5;
Pf = 0.01;
snr_db = 6:2:10; % Signal-to-noise ratio (SNR) in dB
snr = 10.^(snr_db./10); % SNR in linear scale
thresh_theory = 2*gammaincinv(1-Pf,L/2);
u = L./2; % Time-Bandwidth product
PD = [];

for ii = 1:length(snr_db)
    temp1 = (L/2)*snr(ii);
    A = temp1./(1 + temp1);
    n = 0:1:u-2;
    t1 = sum((1./factorial(n)).*(thresh_theory./2).^(n));
    t2 = sum((1./factorial(n)).*(((thresh_theory./2).*(A)).^(n)));
    Pd_theory = exp(-thresh_theory./2).*t1 + (1./A).^(u-1).*(exp(-thresh_theory./(2.*(1+temp1))) - exp(-thresh_theory./2).*t2);
    pm_theory = 1-Pd_theory;
    PD = 0.5 + (Pd_theory-0.5)*sqrt(snr(ii)/(2+snr(ii)));
    temp = 1:1:60;

    for jj=1:1:60
        sum1 = 0;
        sum2 = 0;
        for kk = 1:1:jj
           sum1 = sum1 + PD;
           sum2 = sum2 + (PD*(1-PD));          
        end
        qm(jj) = 1-qfunc((sum1-(jj))/(sqrt(sum2)));
    end
    hold on;
    grid on;
    plot(temp,qm,'-*');
end

legend('SNR = 10dB','SNR = 8dB','SNR = 6dB');
hold off
xlabel('Number of SVUs (N)');
ylabel('Gloal Probability of Misdetection (Qm)');
title('Qm VS Number of SVUs');