%Comparision of Pm vs SNR for time variant and time invariant channel 

clc
close all
clear all

L = 50;
iter =10^5;
Pf = 0.01;
snr_db = -10:1:10;
snr = 10.^(snr_db./10);
%snr = snrLinear*L;
threshold = 2*gammaincinv(1-Pf,L/2);
n=randn(iter,L);

%For time 
for tt = 1:length(snr_db)
    s = []; % Initializtion
    h = []; % Initializtion
    e1 = [];
    mes = randi([0 1],iter,L);
    s = (2.*(mes)-1);
    h1 = (randn(iter,1)+1j*randn(iter,1))./(sqrt(2));
    h = repmat(h1,1,L);
    y = sqrt(snr(tt)).*abs(h).*s + n;
    energy = (abs(y).^2);
    for kk=1:iter
        e1(kk) =sum(energy(kk,:));
    end
    temp1 = (e1 >= threshold);
    i = sum(temp1); 
    Pd_sim(tt) = i/kk; 
    Pm_sim(tt) = 1-Pd_sim(tt);
end

channel = Rayleigh_fading(100,L,150,3,0.00000001);
channel = channel';
h_temp = repmat(channel,100000,1);
%h_temp = h_temp';

%For Time Varying Channel
for tt = 1:length(snr_db)
    mes = [];
    s = []; 
    e2 = [];
    y = [];
    energy = [];
    temp2 = [];
    i=0;
    mes = randi([0 1],iter,L); 
    s = (2.*(mes)-1); 
    y = sqrt(snr(tt)).*abs(h_temp).*s + n;
    energy = (abs(y).^2);
    for kk=1:iter 
        e2(kk) = sum(energy(kk,:)); 
    end
    temp2 = (e2 >= threshold); 
    i = sum(temp2); 
    Pd1(tt) = i/kk;
    Pm1(tt) = 1-Pd1(tt);
end

plot(snr_db,Pm_sim,'-o','LineWidth',2,'color','red')
hold on
grid on
plot(snr_db,Pm1,'-o','LineWidth',2,'color','green')
legend('Time Invarying Channel','Time Varying Channel');
xlabel('SNR (dBs)')
ylabel('Probability Of Miss Detection (Pm)')
title('Pm vs SNR')