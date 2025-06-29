function BSR_score=SE_EE
clc;clear all;close all
rng(0)
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

N_ttd0=ceil(Nt/2*B/fc);
N_ttd1=N_ttd0;
while mod(Nt,N_ttd1)
    N_ttd1=N_ttd1+1;
end
episilon=1e-4;

Pb_dBm=30;
Pb=10^(Pb_dBm/10);%mW

SNRx_dB=20;
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

delta=0.1;
params_pga.acc=1e-4;
params_pga.iter_max=300;%parameters for pga

Iter_max=200;
params_ts.iter_max=Iter_max;
params_ts.local_iter_max=10;%TS stop criteria
params_ts.Len_tl=Iter_max;% length of tabu list

params_TbRs.iter_max=Iter_max;
params_TbRs.local_iter_max=10;
params_TbRs.Len_tl=Iter_max;


BSR_score=AtnDis*Nt*Frac_Bw/8;
avr=100;

%% 
IniMatrix=zeros(avr,1);
SE_DBF=IniMatrix;
SE_DBFtx=IniMatrix;

SE_FTTD16=IniMatrix;
SE_FTTD64=IniMatrix;
SE_FTTD256=IniMatrix;

SE_TTD2_4b=IniMatrix;
SE_TTD2_2b=IniMatrix;
SE_TTD2_1b=IniMatrix;

SE_TTD4_4b=IniMatrix;
SE_TTD4_2b=IniMatrix;
SE_TTD4_1b=IniMatrix;

SE_TTD16_4b=IniMatrix;
SE_TTD16_2b=IniMatrix;
SE_TTD16_1b=IniMatrix;

SE_PS_LSAA_4b=IniMatrix;
SE_PS_LSAA_2b=IniMatrix;
SE_PS_LSAA_1b=IniMatrix;

SW_Apga_TS=IniMatrix;  
Flag_Identity_Ini=1;  
obj_fun_pg=@compute_rate_v0;
obj_fun=@compute_rate_v0;        
parfor n=1:avr
    
    [Hcc,ch_params]=Channel_Gen_compact(Nt,Nr,L,AtnDis,B, ofdm_paras);

        [SE_DBF(n),~,DBF,Proj_Comb]=DBF_water_filling(Hcc,Pb,Ns,noise_power);

         [SE_PS_LSAA_1b(n),~]=HBF_LSAA(Hcc,Pb,Ns,N_rf,noise_power,episilon,PhaseSet_1bit); 
         [SE_PS_LSAA_2b(n),~]=HBF_LSAA(Hcc,Pb,Ns,N_rf,noise_power,episilon,PhaseSet_2bit);   
         [SE_PS_LSAA_4b(n),~]=HBF_LSAA(Hcc,Pb,Ns,N_rf,noise_power,episilon,PhaseSet_4bit); 

         Ndt=1;Ndr=Ndt;
         [SE_TTD2_1b(n),~]=HBF_DDP(Hcc,Pb,Ns,N_rf,Ndt,Ndr,noise_power,AtnDis,B,ch_params,ofdm_paras,PhaseSet_1bit);
         [SE_TTD2_2b(n),~]=HBF_DDP(Hcc,Pb,Ns,N_rf,Ndt,Ndr,noise_power,AtnDis,B,ch_params,ofdm_paras,PhaseSet_2bit);
         [SE_TTD2_4b(n),~]=HBF_DDP(Hcc,Pb,Ns,N_rf,Ndt,Ndr,noise_power,AtnDis,B,ch_params,ofdm_paras,PhaseSet_4bit);
         Ndt=4;Ndr=Ndt;
         [SE_TTD4_1b(n),~]=HBF_DDP(Hcc,Pb,Ns,N_rf,Ndt,Ndr,noise_power,AtnDis,B,ch_params,ofdm_paras,PhaseSet_1bit);
         [SE_TTD4_2b(n),~]=HBF_DDP(Hcc,Pb,Ns,N_rf,Ndt,Ndr,noise_power,AtnDis,B,ch_params,ofdm_paras,PhaseSet_2bit);
         [SE_TTD4_4b(n),~]=HBF_DDP(Hcc,Pb,Ns,N_rf,Ndt,Ndr,noise_power,AtnDis,B,ch_params,ofdm_paras,PhaseSet_4bit)
         Ndt=16;Ndr=Ndt;
         [SE_TTD16_1b(n),~]=HBF_DDP(Hcc,Pb,Ns,N_rf,Ndt,Ndr,noise_power,AtnDis,B,ch_params,ofdm_paras,PhaseSet_1bit);
         [SE_TTD16_2b(n),~]=HBF_DDP(Hcc,Pb,Ns,N_rf,Ndt,Ndr,noise_power,AtnDis,B,ch_params,ofdm_paras,PhaseSet_2bit);
         [SE_TTD16_4b(n),~]=HBF_DDP(Hcc,Pb,Ns,N_rf,Ndt,Ndr,noise_power,AtnDis,B,ch_params,ofdm_paras,PhaseSet_4bit);

         Nfdt=2;Nfdr=Nfdt;
         [SE_FTTD32(n),~]=HBF_FTTD(Hcc,DBF,N_rf,Nfdt,Nfdr,noise_power,AtnDis,B,ofdm_paras);
         Nfdt=8;Nfdr=Nfdt;
         [SE_FTTD64(n),~]=HBF_FTTD(Hcc,DBF,N_rf,Nfdt,Nfdr,noise_power,AtnDis,B,ofdm_paras);
         Nfdt=32;Nfdr=Nfdt;
         [SE_FTTD128(n),~]=HBF_FTTD(Hcc,DBF,N_rf,Nfdt,Nfdr,noise_power,AtnDis,B,ofdm_paras);
               
         [SW_Apga_TS(n),~]=HBF_APGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,noise_power,episilon,delta,obj_fun,obj_fun_pg,Flag_Identity_Ini,Proj_Comb);
   
end
   
time_mark=datestr(now,'mmmm-dd');
file_name=['SE_EE','_',time_mark];
save(file_name,'-v7.3')
len_avr=avr;
% x_num=7;
SE_DBF_avr=mean(SE_DBF(1:len_avr));
SE_FTTD32_avr=mean(SE_FTTD32(1:len_avr));
SE_FTTD64_avr=mean(SE_FTTD64(1:len_avr));
SE_FTTD128_avr=mean(SE_FTTD128(1:len_avr));

SE_PS_LSAA_1b_avr=mean(SE_PS_LSAA_1b(1:len_avr));
SE_PS_LSAA_2b_avr=mean(SE_PS_LSAA_2b(1:len_avr));
SE_PS_LSAA_4b_avr=mean(SE_PS_LSAA_4b(1:len_avr));

SE_TTD2_1b_avr=mean(SE_TTD2_1b(1:len_avr));
SE_TTD2_2b_avr=mean(SE_TTD2_2b(1:len_avr));
SE_TTD2_4b_avr=mean(SE_TTD2_4b(1:len_avr));

SE_TTD4_1b_avr=mean(SE_TTD4_1b(1:len_avr));
SE_TTD4_2b_avr=mean(SE_TTD4_2b(1:len_avr));
SE_TTD4_4b_avr=mean(SE_TTD4_4b(1:len_avr));

SE_TTD16_1b_avr=mean(SE_TTD16_1b(1:len_avr));
SE_TTD16_2b_avr=mean(SE_TTD16_2b(1:len_avr));
SE_TTD16_4b_avr=mean(SE_TTD16_4b(1:len_avr));

SW_Apga_TS_avr=mean(SW_Apga_TS(1:len_avr)); 


EE_DBF=EE_compute_DC(SE_DBF_avr,Nr,Nt);

EE_FTTD32=EE_compute_FTTD_HBF(SE_FTTD32_avr,Nr,Nt,N_rf,2,2);
EE_FTTD64=EE_compute_FTTD_HBF(SE_FTTD64_avr,Nr,Nt,N_rf,8,8);
EE_FTTD128=EE_compute_FTTD_HBF(SE_FTTD128_avr,Nr,Nt,N_rf,32,32);

EE_PS_LSAA_1b=EE_compute_PS_HBF(SE_PS_LSAA_1b_avr,Nr,Nt,N_rf,1);
EE_PS_LSAA_2b=EE_compute_PS_HBF(SE_PS_LSAA_2b_avr,Nr,Nt,N_rf,2);
EE_PS_LSAA_4b=EE_compute_PS_HBF(SE_PS_LSAA_4b_avr,Nr,Nt,N_rf,4);

EE_TTD2_1b=EE_compute_TTD_HBF(SE_TTD2_1b_avr,Nr,Nt,N_rf,1,1,1);
EE_TTD2_2b=EE_compute_TTD_HBF(SE_TTD2_2b_avr,Nr,Nt,N_rf,1,1,2);
EE_TTD2_4b=EE_compute_TTD_HBF(SE_TTD2_4b_avr,Nr,Nt,N_rf,1,1,4);

EE_TTD4_1b=EE_compute_TTD_HBF(SE_TTD4_1b_avr,Nr,Nt,N_rf,4,4,1);
EE_TTD4_2b=EE_compute_TTD_HBF(SE_TTD4_2b_avr,Nr,Nt,N_rf,4,4,2);
EE_TTD4_4b=EE_compute_TTD_HBF(SE_TTD4_4b_avr,Nr,Nt,N_rf,4,4,4);

EE_TTD16_1b=EE_compute_TTD_HBF(SE_TTD16_1b_avr,Nr,Nt,N_rf,16,16,1);
EE_TTD16_2b=EE_compute_TTD_HBF(SE_TTD16_2b_avr,Nr,Nt,N_rf,16,16,2);
EE_TTD16_4b=EE_compute_TTD_HBF(SE_TTD16_4b_avr,Nr,Nt,N_rf,16,16,4);


SW_apga_EE=EE_compute_SW_HBF(SW_Apga_TS_avr,Nr,Nt,N_rf);


time_mark=datestr(now,'mmmm-dd');
file_name=['SE_EE','_',time_mark];
save(file_name,'-v7.3')
   
altstyles = {'m--o'; 'r-p';'g-h'; 'b-s'; 'c-v'; 'k-*'; 'r-.p';'g-.h'; 'b-.s'; 'c-.v'; 'k-.*'};


figure
hold on
h0=plot(SE_DBF_avr,EE_DBF,'k*','LineWidth',1.5)
h1=plot(SW_Apga_TS_avr,SW_apga_EE,'rv','LineWidth',1.5);
h2=plot(SE_PS_LSAA_1b_avr,EE_PS_LSAA_1b,'p','LineWidth',1.5);
h3=plot(SE_PS_LSAA_2b_avr,EE_PS_LSAA_2b,'s','LineWidth',1.5);
h4=plot(SE_PS_LSAA_4b_avr,EE_PS_LSAA_4b,'x','LineWidth',1.5);

h5=plot(SE_FTTD32_avr,EE_FTTD32,'v','LineWidth',1.5);
h6=plot(SE_FTTD64_avr,EE_FTTD64,'v','LineWidth',1.5);
h7=plot(SE_FTTD128_avr,EE_FTTD128,'v','LineWidth',1.5);

h8=plot(SE_TTD2_1b_avr,EE_TTD2_1b,'bp','LineWidth',1.5);
h9=plot(SE_TTD2_2b_avr,EE_TTD2_2b,'bs','LineWidth',1.5);
h10=plot(SE_TTD2_4b_avr,EE_TTD2_4b,'bx','LineWidth',1.5);


h11=plot(SE_TTD4_1b_avr,EE_TTD4_1b,'bp','LineWidth',1.5);
h12=plot(SE_TTD4_2b_avr,EE_TTD4_2b,'bs','LineWidth',1.5);
h13=plot(SE_TTD4_4b_avr,EE_TTD4_4b,'bx','LineWidth',1.5);

h14=plot(SE_TTD16_1b_avr,EE_TTD16_1b,'bp','LineWidth',1.5);
h15=plot(SE_TTD16_2b_avr,EE_TTD16_2b,'bs','LineWidth',1.5);
h16=plot(SE_TTD16_4b_avr,EE_TTD16_4b,'bx','LineWidth',1.5);

legend([h0,h1,h2,h3,h4],'DBF','SW-HBF','PS-HBF:1bits','PS-HBF:2bits','PS-HBF:4bits' );
xlabel('Average SE (bits/s/Hz)')
ylabel('Average EE (Gbits/J)')
grid on
box on

axesNew = axes('position',get(gca,'position'),'visible','off');
legend(axesNew,[h5,h8,h9,h10],'FTTD-HBF','TTD-HBF:1bits','TTD-HBF:2bits','TTD-HBF:4bit','location','west')

set(h2,'Color',[0.6,0.2,0])%brown
set(h3,'Color',[0.6,0.2,0])%brown
set(h4,'Color',[0.6,0.2,0])%brown

set(h5,'Color',[0,0.5,0])%牛油绿
set(h6,'Color',[0,0.5,0])%牛油绿
set(h7,'Color',[0,0.5,0])%牛油绿
% figure
% hold on
% h0=plot(SE_DBF_avr,EE_DBF,'k*','LineWidth',1.5)
% h1=plot(SW_Apga_TS_avr,SW_apga_EE,'bv','LineWidth',1.5);
% h2=plot(SE_PS_LSAA_1b_avr,EE_PS_LSAA_1b,'mp','LineWidth',1.5);
% h3=plot(SE_PS_LSAA_2b_avr,EE_PS_LSAA_2b,'ms','LineWidth',1.5);
% h4=plot(SE_PS_LSAA_4b_avr,EE_PS_LSAA_4b,'mx','LineWidth',1.5);
% 
% h5=plot(SE_FTTD32_avr,EE_FTTD32,'gh','LineWidth',1.5);
% h6=plot(SE_FTTD64_avr,EE_FTTD64,'gh','LineWidth',1.5);
% h7=plot(SE_FTTD128_avr,EE_FTTD128,'gh','LineWidth',1.5);
% 
% h8=plot(SE_TTD2_1b_avr,EE_TTD2_1b,'rp','LineWidth',1.5);
% h9=plot(SE_TTD2_2b_avr,EE_TTD2_2b,'rs','LineWidth',1.5);
% h10=plot(SE_TTD2_4b_avr,EE_TTD2_4b,'rx','LineWidth',1.5);
% 
% 
% h11=plot(SE_TTD4_1b_avr,EE_TTD4_1b,'rp','LineWidth',1.5);
% h12=plot(SE_TTD4_2b_avr,EE_TTD4_2b,'rs','LineWidth',1.5);
% h13=plot(SE_TTD4_4b_avr,EE_TTD4_4b,'rx','LineWidth',1.5);
% 
% h14=plot(SE_TTD16_1b_avr,EE_TTD16_1b,'rp','LineWidth',1.5);
% h15=plot(SE_TTD16_2b_avr,EE_TTD16_2b,'rs','LineWidth',1.5);
% h16=plot(SE_TTD16_4b_avr,EE_TTD16_4b,'rx','LineWidth',1.5);
% 
% legend([h0,h1,h2,h3,h4,h5,h8,h9,h10],'DBF','SW-HBF','PS-HBF:1bits','PS-HBF:2bits','PS-HBF:4bits',...
%     'FTTD-HBF','TTD-HBF:1bits','TTD-HBF:2bits','TTD-HBF:4bit');
% xlabel('Average SE (bits/s/Hz)')
% ylabel('Average EE (Gbits/J)')
% grid on
% box on

ccc=1;







