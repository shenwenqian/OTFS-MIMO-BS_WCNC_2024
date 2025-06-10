clear
clc
%% Parameter settings
M  = 4;
N  = 4;% number of DD unit
fc = 3*10e9;%carrier Frequency
fs = 15*10e6;%carrier interval
Ntx= 8;%number of transmit antenna
Nrx= 8;%number of receive antenna
L  = 4; %number of path
c  = 3*10e8;%speed of light 
lamada = c/fc; % wavelength
d = lamada/2; % 信道间隔

SNR_dB = -10:5:15;
SNR_table = 10.^(SNR_dB/10);%SNR
Ns = 4; %data streams 
Ptot = 10; %Power
count = 100;

for i_SNR = 1:length(SNR_table)

    i_SNR
    SNR = SNR_table(i_SNR);
    Q1 = 0;
    Q2 = 0;
    Q3 = 0;
    for i_count=1:count

        [Hbs,Hqt] = channel(M,N,L,Nrx,Ntx,c,fc,fs,d); %channel generation 

        Q1 = digitalEX(Hqt,SNR,Ntx,Ns,M,N,Ptot)+Q1;
        Q2 = digitalEXerror(Hqt,Hbs,SNR,Ntx,Ns,M,N,Ptot)+Q2;
        Q3 = digitalone1error(Hqt,Hbs,SNR,Ntx,Nrx,M,N,Ptot)+Q3;

    end

    C1(i_SNR)=real(10*Q1/M/N/Ns/count);
    C2(i_SNR)=real(10*Q2/M/N/Ns/count);
    C3(i_SNR)=real(10*Q3/M/N/Ns/count);
end

plot(SNR_dB,C1,'-or','Linewidth',1.5);
hold on;
plot(SNR_dB,C2,'--og','Linewidth',1.5);
hold on;
plot(SNR_dB,C3,'--*b','Linewidth',1.5);
hold on;


legend('Proposed scheme with BSC','Proposed scheme without BSC','Traditional scheme');
xlabel('SNR (dB)');
ylabel('Achievable rate (bps/Hz)');

           
    
        











