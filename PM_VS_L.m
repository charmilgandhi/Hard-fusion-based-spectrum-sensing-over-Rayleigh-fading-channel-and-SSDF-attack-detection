%Pm vs Number of Samples (L)

clear all;
close all;
clc;
L = 0:10:500;
L(1) = 1;
pf = 0.1;
threshold = qfuncinv(pf);
for snr = -16:2:-4
    snrL = db2pow(snr);
    for i = 1:length(L)
        snrLinear = snrL*L(i);
        Pmi(i) = qfunc((snrLinear - (sqrt(L(i))*threshold))/(sqrt(L(i)+2*snrLinear)));
        %PMi(i) = 1/2 + (Pmi(i) - 1/2)*(sqrt((snrLinear)/(2+snrLinear)));
        A = 1/(sqrt(L(i)+2*snrLinear));
        B = (sqrt(L(i))*threshold)/(sqrt(L(i)+2*snrLinear));
        PMiBar(i) = (qfunc(-1*B) - ((exp(-B/(A*snrLinear)))*(exp(1/(2*(A.*A)*(snrLinear*snrLinear))))*(qfunc((1/(A*snrLinear)) - B))));
    end
    plot(L, Pmi,'-x','LineWidth',1.5);
    hold on
end

grid on
legend('-16dB', '-14dB','-12dB','-10dB', '-8dB', '-6dB', '-4dB');
axis([0 500 0 1]);
xlabel('Number Of Samples (L)')
ylabel('Probability Of Miss Detection (For ith SVU) (Pm)')
title('Pm vs L for different values of SNR');