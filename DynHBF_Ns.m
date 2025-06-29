function TC=DybHBF_Ns
% clc;clear all;close all
tic
rng(0)
GHz=1e+9;%bandwidth=1GHz
Nc=64;% carrier number
D=Nc/4;% maximum delay
B=30*GHz;% bandwidth
fc=300*GHz;% carrier frequency
Nt=256;
Nr=256;
% N_rf=5;
% Ns=N_rf;

AtnDis=1/2;
Frac_Bw=B/fc;

ofdm_paras.B=B;
ofdm_paras.Nc=Nc;% carrier number
ofdm_paras.D=D;% maximum delay
ofdm_paras.fc=fc;% carrier frequency
L=4;

episilon=1e-4;

Pb_dBm=30;
Pb=10^(Pb_dBm/10);%mW single carrier

SNRx_dB=20;
% SNRx_dB=-10:5:30;
SNRx=10.^(SNRx_dB/10);
noise_power=Pb./SNRx;% noise power

Phase_Bit=[1 2 4];
phase_bit_num=length(Phase_Bit);
Nq_set=2.^Phase_Bit;
for lb=1:phase_bit_num   
phase_Set{lb}=(-pi:2*pi/Nq_set(lb):pi);
end
PhaseSet_1bit=phase_Set{1};
PhaseSet_2bit=phase_Set{2};
PhaseSet_4bit=phase_Set{3};

Ns_all=1:5;
LenN=length(Ns_all);

BSR_score=AtnDis*Nt*Frac_Bw/8;
MC=100;

iter_max=500;tol_obj=1e-4;
%% 
IniMatrix=zeros(MC,LenN);
SE_DBF=IniMatrix;
SE_PS_LSAA_4b=IniMatrix;
SE_PS_Dyn_4b=IniMatrix;

SE_PS_LSAA_2b=IniMatrix;
SE_PS_Dyn_2b=IniMatrix;

SE_PS_LSAA_1b=IniMatrix;
SE_PS_Dyn_1b=IniMatrix;

parfor ii=1:MC
    [Hcc,ch_params]=Channel_Gen_compact(Nt,Nr,L,AtnDis,B,ofdm_paras);
    for jj=1:LenN
        Ns=Ns_all(jj);
        N_rf=Ns;


%         [SE_DBF(ii,jj),~,DBF,Proj_Comb]=DBF_water_filling(Hcc,Pb,Ns,noise_power);
% 
%         [SE_PS_LSAA_4b(ii,jj),~]=HBF_pp_quantized(Hcc,Pb,Ns,N_rf,noise_power,episilon,PhaseSet_4bit);
%         [SE_PS_LSAA_2b(ii,jj),~]=HBF_pp_quantized(Hcc,Pb,Ns,N_rf,noise_power,episilon,PhaseSet_2bit);
%         [SE_PS_LSAA_1b(ii,jj),~]=HBF_pp_quantized(Hcc,Pb,Ns,N_rf,noise_power,episilon,PhaseSet_1bit);

        SE_PS_Dyn_4b(ii,jj)=Dyn_PS_HBF(Hcc,N_rf,Ns,Pb,noise_power,PhaseSet_4bit,iter_max,tol_obj);        
        SE_PS_Dyn_2b(ii,jj)=Dyn_PS_HBF(Hcc,N_rf,Ns,Pb,noise_power,PhaseSet_2bit,iter_max,tol_obj);              
        SE_PS_Dyn_1b(ii,jj)=Dyn_PS_HBF(Hcc,N_rf,Ns,Pb,noise_power,PhaseSet_1bit,iter_max,tol_obj);

    end
end
TC=toc;


SE_DBF_avr=mean(SE_DBF);
SE_PS_LSAA_1b_avr=mean(SE_PS_LSAA_1b);
SE_PS_LSAA_2b_avr=mean(SE_PS_LSAA_2b);
SE_PS_LSAA_4b_avr=mean(SE_PS_LSAA_4b);

SE_PS_Dyn_1b_avr=mean(SE_PS_Dyn_1b);
SE_PS_Dyn_2b_avr=mean(SE_PS_Dyn_2b);
SE_PS_Dyn_4b_avr=mean(SE_PS_Dyn_4b);

time_mark=datestr(now,'mmmm-dd');
file_name=['DynHBF_Ns_',time_mark];
save(file_name,'-v7.3')
% 
% altstyles = {'m--o'; 'r-p';'g-h'; 'b-s'; 'c-v'; 'k-*'; 'r-.p';'g-.h'; 'b-.s'; 'c-.v'; 'k-.*'};
% % x=B_all/GHz;
% xx=SNRx_dB;
% figure
% hold on
% h0=plot(xx,SE_DBF_avr,'r-*','LineWidth',1.5);
% h1=plot(xx,SE_PS_LSAA_4b_avr,'b-o','LineWidth',1.5);
% h2=plot(xx,SE_PS_LSAA_2b_avr,'b--o','LineWidth',1.5);
% h3=plot(xx,SE_PS_LSAA_1b_avr,'b:o','LineWidth',1.5);
% h4=plot(xx,SE_PS_Dyn_4b_avr,'k-s','LineWidth',1.5);
% h5=plot(xx,SE_PS_Dyn_2b_avr,'k--s','LineWidth',1.5);
% h6=plot(xx,SE_PS_Dyn_1b_avr,'k:s','LineWidth',1.5);
% 
% % set(h1,'Color',[0,0.5,0])%牛油绿
% % set(h2,'Color',[0.6,0.2,0])%brown
% % set(h3,'Color',[0.6,0.2,0])%brown
% % set(h4,'Color',[0.6,0.2,0])%brown
% 
% % set(h11,'Color',[0.4940 0.1840 0.5560])% purple
% 
% legend([h0,h1,h2,h3,h4,h5,h6],'DBF','FC-PS-HBF:4bits','FC-PS-HBF:2bits','FC-PS-HBF:1bit', ...
%     'Dyn-PS-HBF:4bits','Dyn-PS-HBF:2bits','Dyn-PS-HBF:1bit');
% ylabel('Average SE (bits/s/Hz)')
% xlabel('SNR(dB)')
% % xticks(Nall)
% xlim([xx(1) xx(end)])
% grid on
% box on



ccc=1;
