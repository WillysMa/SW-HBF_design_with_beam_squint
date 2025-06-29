clc;clear all;close all
%% variable: SNR
GHz=1e+9;%bandwidth=1GHz
Nc=64;% carrier number
D=Nc/4;% maximum delay
B=30*GHz;% bandwidth
fc=300*GHz;% carrier frequency
Nt=256;
Nr=Nt;
N_rf=4;
Ns=4;
AtnDis=1/2;
Frac_Bw=B/fc;

ofdm_paras.B=B;
ofdm_paras.Nc=Nc;% carrier number
ofdm_paras.D=D;% maximum delay
ofdm_paras.fc=fc;% carrier frequency
L=4;


Pb_dBm=30;
Pb=10^(Pb_dBm/10);%mW

SNRx_dB=10;
SNRx=10.^(SNRx_dB/10);
noise_power=Pb./SNRx;% noise power

BSR_score=AtnDis*Nt*Frac_Bw/8;

Nfttd_all=2.^(2:8);
len=length(Nfttd_all);
avr=2;
SE_all=zeros(avr,len);
SE_all_tx=zeros(avr,len);
for n=1:avr
    
    [Hcc,ch_params]=Channel_Gen_compact(Nt,Nr,L,AtnDis,B, ofdm_paras);

        [SE_DBF(n),~,DBF]=DBF_water_filling(Hcc,Pb,Ns,noise_power);

    for id=1:len
         Nfdt=Nfttd_all(id);Nfdr=Nfdt;
         [SE_all(n,id),SE_all_tx(n,id)]=HBF_FTTD(Hcc,DBF,N_rf,Nfdt,Nfdr,noise_power,AtnDis,B,ofdm_paras);
    end
   
end

SE_FTTD=mean(SE_all);
SE_FTTD_tx=mean(SE_all_tx);
time_mark=datestr(now,'mmmm-dd');
file_name=['FTTD_test','_',time_mark];
save(file_name,'-v7.3')

figure
plot(Nfttd_all,SE_FTTD_tx,'b-o','LineWidth',1.5)
xlabel('Number of FTTD each RF chain')
ylabel('Average SE (bits/s/Hz)')
grid on

figure
plot(Nfttd_all,SE_FTTD,'b-o','LineWidth',1.5)
xlabel('Number of FTTD each RF chain')
ylabel('Average SE (bits/s/Hz)')
grid on
ccc=1;
