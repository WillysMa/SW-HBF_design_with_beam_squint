function TC=SE_SNR(Ndt)
clc;clear all;close all
%% variable: SNR
rng(0)
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
% while mod(Nt,N_ttd1)
%     N_ttd1=N_ttd1+1;
% end
if (~exist('Ndt','var'))
    Ndt=1;
end
Ndr=Ndt;
episilon=1e-4;

% Nfdt=Nt/N_rf;
% Nfdr=Nr/N_rf;
Nfdt=2;
Nfdr=2;

saved_name=['SE_SNRttd',num2str(Ndt)];

Pb_dBm=30;
Pb=10^(Pb_dBm/10);%mW single carrier

SNRx_dB=20:5:30;
SNRx=10.^(SNRx_dB/10);
noise_power_all=Pb./SNRx;% noise power


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
params_pga.iter_max=1000;%parameters for pga

Iter_max=200;
params_ts.iter_max=Iter_max;
params_ts.local_iter_max=10;%TS stop criteria
params_ts.Len_tl=Iter_max;% length of tabu list

params_TbRs.iter_max=Iter_max;
params_TbRs.local_iter_max=10;
params_TbRs.Len_tl=Iter_max;

len=length(SNRx_dB);

BSR_score=AtnDis*Nt*Frac_Bw/8;
avr=100;
TC=zeros(1,avr);
%% 
IniMatrix=zeros(avr,len);
SE_DBF=IniMatrix;
SE_DBFtx=IniMatrix;

SE_FTTD=IniMatrix;
SE_FTTD0=IniMatrix;
SE_TTD_4b=IniMatrix;
SE_TTD_4b_tx=IniMatrix;

SE_TTD_2b=IniMatrix;
SE_TTD_1b=IniMatrix;

SE_PS_LSAA_4b=IniMatrix;
SE_PS_LSAA_4b_tx=IniMatrix;

SE_PS_LSAA_2b=IniMatrix;
SE_PS_LSAA_1b=IniMatrix;

SW_rs=IniMatrix;
SW_Rdn_TS=IniMatrix; 
SW_Apga_TS=IniMatrix;

SW_rs1=IniMatrix;
SW_Rdn_TS1=IniMatrix; 
SW_Apga_TS1=IniMatrix;


SW_random=IniMatrix; 

obj_fun_pg=@compute_rate_v0;
obj_fun=@compute_rate_v0;     
tic
for n=1:avr
[Hcc,ch_params]=Channel_Gen_compact(Nt,Nr,L,AtnDis,B,ofdm_paras);

%     tic
    for id=1:len   
        noise_power=noise_power_all(id);
        [SE_DBF(n,id),~,DBF,Proj_Comb]=DBF_water_filling(Hcc,Pb,Ns,noise_power);


         [SE_PS_LSAA_1b(n,id),~]=HBF_pp_quantized(Hcc,Pb,Ns,N_rf,noise_power,episilon,PhaseSet_1bit); 
         [SE_PS_LSAA_2b(n,id),~]=HBF_pp_quantized(Hcc,Pb,Ns,N_rf,noise_power,episilon,PhaseSet_2bit);   
         [SE_PS_LSAA_4b(n,id),~]=HBF_pp_quantized(Hcc,Pb,Ns,N_rf,noise_power,episilon,PhaseSet_4bit); 
   
         [SE_PS_LSAA_1b(n,id),~]=HBF_LSAA(Hcc,Pb,Ns,N_rf,noise_power,episilon,PhaseSet_1bit); 
%          [SE_PS_LSAA_2b(n,id),~]=HBF_LSAA(Hcc,Pb,Ns,N_rf,noise_power,episilon,PhaseSet_2bit);   
%          [SE_PS_LSAA_4b(n,id),~]=HBF_LSAA(Hcc,Pb,Ns,N_rf,noise_power,episilon,PhaseSet_4bit); 

         [SE_TTD_1b(n,id),~]=HBF_DDP(Hcc,Pb,Ns,N_rf,Ndt,Ndr,noise_power,AtnDis,B,ch_params,ofdm_paras,PhaseSet_1bit);
         [SE_TTD_2b(n,id),~]=HBF_DDP(Hcc,Pb,Ns,N_rf,Ndt,Ndr,noise_power,AtnDis,B,ch_params,ofdm_paras,PhaseSet_2bit);
         [SE_TTD_4b(n,id),~]=HBF_DDP(Hcc,Pb,Ns,N_rf,Ndt,Ndr,noise_power,AtnDis,B,ch_params,ofdm_paras,PhaseSet_4bit);

         [SE_FTTD(n,id),~]=HBF_FTTD(Hcc,DBF,N_rf,Nfdt,Nfdr,noise_power,AtnDis,B,ofdm_paras);

         Flag_Identity_Ini=1;         
         [SW_Apga_TS(n,id),~]=HBF_APGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,noise_power,episilon,delta,obj_fun,obj_fun_pg,Flag_Identity_Ini,Proj_Comb);
%          [SW_Rdn_TS(n,id),~]=HBF_TS(params_ts,Pb,Ns,N_rf,Hcc,noise_power,episilon,obj_fun,Flag_Identity_Ini,Proj_Comb);      
         [SW_rs(n,id),~]=HBF_APGA_TbRs(params_pga,params_TbRs,Pb,Ns,N_rf,Hcc,noise_power,episilon,delta,obj_fun,obj_fun_pg,Flag_Identity_Ini,Proj_Comb); 

%          Flag_Identity_Ini=0;         
%          [SW_Apga_TS1(n,id),~]=HBF_APGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,noise_power,episilon,delta,obj_fun,obj_fun_pg,Flag_Identity_Ini,Proj_Comb);
%          [SW_Rdn_TS1(n,id),~]=HBF_TS(params_ts,Pb,Ns,N_rf,Hcc,noise_power,episilon,obj_fun,Flag_Identity_Ini,Proj_Comb);      
%          [SW_rs1(n,id),~]=HBF_APGA_TbRs(params_pga,params_TbRs,Pb,Ns,N_rf,Hcc,noise_power,episilon,delta,obj_fun,obj_fun_pg,Flag_Identity_Ini,Proj_Comb); 

         [SW_random(n,id),~]=SW_HBF_random(Pb,Ns,N_rf,Hcc,noise_power,episilon);         
% 
%       disp(['id=',num2str(id),', avr=',num2str(n)])

        ccc=1;
    end
%     TC(n)=toc/60;
%     disp(['avr=',num2str(n)])
end
toc
time_mark=datestr(now,'mmmm-dd');
file_name=[saved_name,'_',time_mark];
save(file_name,'-v7.3')
len_avr=avr;
% x_num=7;
SE_DBF_avr=mean(SE_DBF(1:len_avr,:));
SE_FTTD_avr=mean(SE_FTTD(1:len_avr,:));
SE_PS_LSAA_1b_avr=mean(SE_PS_LSAA_1b(1:len_avr,:));
SE_PS_LSAA_2b_avr=mean(SE_PS_LSAA_2b(1:len_avr,:));
SE_PS_LSAA_4b_avr=mean(SE_PS_LSAA_4b(1:len_avr,:));

SE_TTD_1b_avr=mean(SE_TTD_1b(1:len_avr,:));
SE_TTD_2b_avr=mean(SE_TTD_2b(1:len_avr,:));
SE_TTD_4b_avr=mean(SE_TTD_4b(1:len_avr,:));

SW_random_avr=mean(SW_random(1:len_avr,:)); 
SW_rs_avr=mean(SW_rs(1:len_avr,:));
SW_Apga_TS_avr=mean(SW_Apga_TS(1:len_avr,:)); 
SW_Rdn_TS_avr=mean(SW_Rdn_TS(1:len_avr,:)); 

% SW_rs_avr1=mean(SW_rs1(1:len_avr,:));
% SW_Apga_TS_avr1=mean(SW_Apga_TS1(1:len_avr,:)); 
% SW_Rdn_TS_avr1=mean(SW_Rdn_TS1(1:len_avr,:)); 

EE_DBF=EE_compute_DC(SE_DBF_avr,Nr,Nt);
EE_FTTD=EE_compute_FTTD_HBF(SE_FTTD_avr,Nr,Nt,N_rf,Nfdt,Nfdr);
EE_PS_LSAA_1b=EE_compute_PS_HBF(SE_PS_LSAA_1b_avr,Nr,Nt,N_rf,1);
EE_PS_LSAA_2b=EE_compute_PS_HBF(SE_PS_LSAA_2b_avr,Nr,Nt,N_rf,2);
EE_PS_LSAA_4b=EE_compute_PS_HBF(SE_PS_LSAA_4b_avr,Nr,Nt,N_rf,4);

EE_TTD_1b=EE_compute_TTD_HBF(SE_TTD_1b_avr,Nr,Nt,N_rf,Ndt,Ndr,1);
EE_TTD_2b=EE_compute_TTD_HBF(SE_TTD_2b_avr,Nr,Nt,N_rf,Ndt,Ndr,2);
EE_TTD_4b=EE_compute_TTD_HBF(SE_TTD_4b_avr,Nr,Nt,N_rf,Ndt,Ndr,4);

EE_SW_rs=EE_compute_SW_HBF(SW_rs_avr,Nr,Nt,N_rf);
EE_SW_random=EE_compute_SW_HBF(SW_random_avr,Nr,Nt,N_rf);
SW_apga_EE=EE_compute_SW_HBF(SW_Apga_TS_avr,Nr,Nt,N_rf);
SW_Rdn_EE=EE_compute_SW_HBF(SW_Rdn_TS_avr,Nr,Nt,N_rf);


time_mark=datestr(now,'mmmm-dd');
%% file_name=['SE_SNR_',time_mark];
file_name=[saved_name,'_',time_mark];
%  save(file_name,'-v7.3')
   
altstyles = {'m--o'; 'r-p';'g-h'; 'b-s'; 'c-v'; 'k-*'; 'r-.p';'g-.h'; 'b-.s'; 'c-.v'; 'k-.*'};
% x=B_all/GHz;
x=SNRx_dB;
figure
hold on
h0=plot(x,SE_DBF_avr,'k-*','LineWidth',1.5);
h1=plot(x,SE_FTTD_avr,'b-v','LineWidth',1.5);
h2=plot(x,SE_PS_LSAA_4b_avr,'-o','LineWidth',1.5);
h3=plot(x,SE_PS_LSAA_2b_avr,'--o','LineWidth',1.5);
h4=plot(x,SE_PS_LSAA_1b_avr,':o','LineWidth',1.5);

h5=plot(x,SE_TTD_4b_avr,'b-s','LineWidth',1.5);
h6=plot(x,SE_TTD_2b_avr,'b--s','LineWidth',1.5);
h7=plot(x,SE_TTD_1b_avr,'b:s','LineWidth',1.5);

h8=plot(x,SW_Apga_TS_avr,'r-d','LineWidth',1.5);
h9=plot(x,SW_Rdn_TS_avr,'r--d','LineWidth',1.5);
h10=plot(x,SW_rs_avr,'r:d','LineWidth',1.5);
h11=plot(x,SW_random_avr,'m-x','LineWidth',1.5);


set(h1,'Color',[0,0.5,0])%牛油绿
set(h2,'Color',[0.6,0.2,0])%brown
set(h3,'Color',[0.6,0.2,0])%brown
set(h4,'Color',[0.6,0.2,0])%brown

% set(h11,'Color',[0.4940 0.1840 0.5560])% purple

legend([h0,h1,h2,h3,h4,h5,h6,h7],'DBF','FTTD-HBF','PS-HBF:4bits','PS-HBF:2bits','PS-HBF:1bit','TTD-HBF:4bits','TTD-HBF:2bits','TTD-HBF:1bit');
ylabel('Average SE (bits/s/Hz)')
xlabel('SNR(dB)')
% xticks(Nall)
xlim([SNRx_dB(1) SNRx_dB(end)])
grid on
box on

axesNew = axes('position',get(gca,'position'),'visible','off');
legend(axesNew,[h8,h9,h10,h11],'SW-HBF:PGA-TS','SW-HBF:TS','SW-HBF:PGA-TbRs','SW-HBF:random','location','west')



figure
hold on
h0=plot(x,EE_DBF,'k-*','LineWidth',1.5);
h1=plot(x,EE_FTTD,'b-v','LineWidth',1.5);
h2=plot(x,EE_PS_LSAA_4b,'m-o','LineWidth',1.5);
h3=plot(x,EE_PS_LSAA_2b,'m--o','LineWidth',1.5);
h4=plot(x,EE_PS_LSAA_1b,'m-.o','LineWidth',1.5);

h5=plot(x,EE_TTD_4b,'-s','LineWidth',1.5);
h6=plot(x,EE_TTD_2b,'--s','LineWidth',1.5);
h7=plot(x,EE_TTD_1b,':s','LineWidth',1.5);

h8=plot(x,SW_apga_EE,'-d','LineWidth',1.5);
h9=plot(x,SW_Rdn_EE,'--d','LineWidth',1.5);
h10=plot(x,EE_SW_rs,':d','LineWidth',1.5);
h11=plot(x,EE_SW_random,'-x','LineWidth',1.5);

set(h5,'Color',[0,0.5,0])%牛油绿
set(h6,'Color',[0,0.5,0])%牛油绿
set(h7,'Color',[0,0.5,0])%牛油绿

set(h8,'Color',[0.6,0.2,0])%brown
set(h9,'Color',[0.6,0.2,0])%brown
set(h10,'Color',[0.6,0.2,0])%brown
set(h11,'Color',[0.4940 0.1840 0.5560])% purple

legend([h0,h1,h2,h3,h4,h5,h6,h7],'DBF','FTTD-HBF','PS-HBF:4bits','PS-HBF:2bits','PS-HBF:1bit','TTD-HBF:4bits','TTD-HBF:2bits','TTD-HBF:1bit');

ylabel('Average EE (Gbits/J)')
xlabel('SNR(dB)')
% xticks(Nall)
xlim([SNRx_dB(1) SNRx_dB(end)])
grid on 
box on
axesNew = axes('position',get(gca,'position'),'visible','off');
legend(axesNew,[h8,h9,h10,h11],'SW-HBF:PGA-TS','SW-HBF:TS','SW-HBF:PGA-TbRs','SW-HBF:random','location','west')

ccc=1;


end



