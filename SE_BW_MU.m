function TC=SE_BW_MU
clc;clear all;close all
%% variable: BW
rng(0)
tic
GHz=1e+9;%bandwidth=1GHz
Nc=64;% carrier number
D=Nc/4;% maximum delay
B=30*GHz;% bandwidth
fc=300*GHz;% carrier frequency
Nt=256;
% Nu=4;
% Nrf=Nu;
AtnDis=1/2;

BSR_coef=AtnDis*Nt/fc/8;
ofdm_paras.B=B;
ofdm_paras.Nc=Nc;% carrier number
ofdm_paras.D=D;% maximum delay
ofdm_paras.fc=fc;% carrier frequency
L=4;

radius=50;
r_h=10;


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

params_pga.acc=1e-4;
params_pga.iter_max=1e+3;%parameters for pga

Iter_max=200;
params_ts.iter_max=Iter_max;
params_ts.local_iter_max=10;%TS stop criteria
params_ts.Len_tl=Iter_max;% length of tabu list

iter_max=1e+3;tol=1e-4;

BSR_all=0.1:0.3:1.8;
% BSR_all=0.1:0.3:0.4;
B_all=BSR_all/BSR_coef;
len=length(B_all);

avr=100;

%% 
IniMatrix=zeros(avr,len);
SE_DBF=IniMatrix;
SE_PS_4b=IniMatrix;
SE_PS_2b=IniMatrix;
SE_PS_1b=IniMatrix;
SW_PgaTS=IniMatrix;  
SW_PgaRdn=IniMatrix;  

SE_DBF_v1=IniMatrix;
SE_PS_4b_v1=IniMatrix;
SE_PS_2b_v1=IniMatrix;
SE_PS_1b_v1=IniMatrix;
SW_PgaTS_v1=IniMatrix;  
SW_PgaRdn_v1=IniMatrix;           
parfor n=1:avr   
    for id=1:len   
        B=B_all(id);
        
        % Hall=channel_MISO(Nt,Nu,L,AtnDis,B,ofdm_paras,large_fading); 

        Nu=4; Nrf=Nu;
        Wgtu=ones(1,Nu);
        Hall=Channel_MU(Nt,Nu,L,AtnDis,B,ofdm_paras);       
        SE_DBF(n,id)=WMMSE_MISO(Hall,Wgtu,Pb,noise_power,50,tol);
        [SW_PgaTS(n,id),~,SW_PgaRdn(n,id)]=TwoStepHBF_design_MISO(params_ts,params_pga,Hall,Wgtu,Nrf,Pb,noise_power,iter_max,tol);
        SE_PS_4b(n,id)=PSHBF_MISO_WeiYu(Hall,Wgtu,Nrf,Pb,noise_power,PhaseSet_4bit,iter_max,tol);
        SE_PS_2b(n,id)=PSHBF_MISO_WeiYu(Hall,Wgtu,Nrf,Pb,noise_power,PhaseSet_2bit,iter_max,tol);
        SE_PS_1b(n,id)=PSHBF_MISO_WeiYu(Hall,Wgtu,Nrf,Pb,noise_power,PhaseSet_1bit,iter_max,tol);        


        Nu=8; Nrf=Nu;
        Wgtu=ones(1,Nu);
        Hall=Channel_MU(Nt,Nu,L,AtnDis,B,ofdm_paras);
        SE_DBF_v1(n,id)=WMMSE_MISO(Hall,Wgtu,Pb,noise_power,50,tol);
        [SW_PgaTS_v1(n,id),~,SW_PgaRdn_v1(n,id)]=TwoStepHBF_design_MISO(params_ts,params_pga,Hall,Wgtu,Nrf,Pb,noise_power,iter_max,tol);
        SE_PS_4b_v1(n,id)=PSHBF_MISO_WeiYu(Hall,Wgtu,Nrf,Pb,noise_power,PhaseSet_4bit,iter_max,tol);
        SE_PS_2b_v1(n,id)=PSHBF_MISO_WeiYu(Hall,Wgtu,Nrf,Pb,noise_power,PhaseSet_2bit,iter_max,tol);
        SE_PS_1b_v1(n,id)=PSHBF_MISO_WeiYu(Hall,Wgtu,Nrf,Pb,noise_power,PhaseSet_1bit,iter_max,tol);  
%       disp(['id=',num2str(id),', avr=',num2str(n)])
        ccc=1;
    end
%     disp(['avr=',num2str(n)])
%     time_mark=datestr(now,'mmmm-dd');
%     file_name=['SE_SNR_MU_',time_mark];
%     save(file_name,'-v7.3')
end
TC=toc;
len_avr=avr;
Nu=4;
Nrf=Nu;
SE_DBF_avr=mean(SE_DBF(1:len_avr,:));
SE_PS_1b_avr=mean(SE_PS_1b(1:len_avr,:));
SE_PS_2b_avr=mean(SE_PS_2b(1:len_avr,:));
SE_PS_4b_avr=mean(SE_PS_4b(1:len_avr,:));
SW_PgaTS_avr=mean(SW_PgaTS(1:len_avr,:));
SW_PgaRdn_avr=mean(SW_PgaRdn(1:len_avr,:)); 

EE_DBF=EE_compute_DC_MU(SE_DBF_avr,Nt);
EE_PS_1b=EE_compute_PS_HBF_MU(SE_PS_1b_avr,Nt,Nrf,1);
EE_PS_2b=EE_compute_PS_HBF_MU(SE_PS_2b_avr,Nt,Nrf,2);
EE_PS_4b=EE_compute_PS_HBF_MU(SE_PS_4b_avr,Nt,Nrf,4);
EE_SW_pgats=EE_compute_SW_HBF_MU(SW_PgaTS_avr,Nt,Nrf);
EE_SW_rdn=EE_compute_SW_HBF_MU(SW_PgaRdn_avr,Nt,Nrf);

Nu=8;
Nrf=Nu;
SE_DBF_avr_v1=mean(SE_DBF_v1(1:len_avr,:));
SE_PS_1b_avr_v1=mean(SE_PS_1b_v1(1:len_avr,:));
SE_PS_2b_avr_v1=mean(SE_PS_2b_v1(1:len_avr,:));
SE_PS_4b_avr_v1=mean(SE_PS_4b_v1(1:len_avr,:));
SW_PgaTS_avr_v1=mean(SW_PgaTS_v1(1:len_avr,:));
SW_PgaRdn_avr_v1=mean(SW_PgaRdn_v1(1:len_avr,:)); 

EE_DBF_v1=EE_compute_DC_MU(SE_DBF_avr_v1,Nt);
EE_PS_1b_v1=EE_compute_PS_HBF_MU(SE_PS_1b_avr_v1,Nt,Nrf,1);
EE_PS_2b_v1=EE_compute_PS_HBF_MU(SE_PS_2b_avr_v1,Nt,Nrf,2);
EE_PS_4b_v1=EE_compute_PS_HBF_MU(SE_PS_4b_avr_v1,Nt,Nrf,4);
EE_SW_pgats_v1=EE_compute_SW_HBF_MU(SW_PgaTS_avr_v1,Nt,Nrf);
EE_SW_rdn_v1=EE_compute_SW_HBF_MU(SW_PgaRdn_avr_v1,Nt,Nrf);
time_mark=datestr(now,'mmmm-dd');
file_name=['SE_BW_MU_',time_mark];
save(file_name,'-v7.3')
   
%% draw figure
lineSpec1={'r-o';'k-s';'b-+';'m-^';'g-d';'c-x';{[0.6,0.2,0],'-*'}};
lineSpec2={'r--o';'k--s';'b--+';'m--^';'g--d';'c--x';{[0.6,0.2,0],'--*'}};
lineSpec3={'r:o';'k:s';'b:+';'m:^';'g:d';'c:x';{[0.6,0.2,0],':*'}};
lineSpec4={'r-.o';'k-.s';'b-.+';'m-.^';'g-.d';'c-.x';{[0.6,0.2,0],'-.*'}};
lineSpec0={'ro';'ks';'b+';'m^';'gd';'cx';{[0.6,0.2,0],'*'}};

xx=B_all/GHz;
% xx=BSR_all;
figure
hold on

h1=plot(xx,SE_DBF_avr,'k-*','LineWidth',1.5);
h12=plot(xx,SE_PS_4b_avr,'--o','LineWidth',1.5);
h2=plot(xx,SE_PS_2b_avr,'-o','LineWidth',1.5);
h3=plot(xx,SE_PS_1b_avr,':o','LineWidth',1.5);
h4=plot(xx,SW_PgaTS_avr,'r-s','LineWidth',1.5);
h5=plot(xx,SW_PgaRdn_avr,'b-^','LineWidth',1.5);

% h1=plot(xx,SE_DBF_avr_v1,'k-*','LineWidth',1.5);
% h12=plot(xx,SE_PS_4b_avr_v1,'-o','LineWidth',1.5);
% h2=plot(xx,SE_PS_2b_avr_v1,':o','LineWidth',1.5);
% h3=plot(xx,SE_PS_1b_avr_v1,'--o','LineWidth',1.5);
% h4=plot(xx,SW_PgaTS_avr_v1,'r-s','LineWidth',1.5);
% h5=plot(xx,SW_PgaRdn_avr_v1,'b-^','LineWidth',1.5);
% set(h5,'Color',[0,0.5,0])%牛油绿
% set(h11,'Color',[0.4940 0.1840 0.5560])% purple
% set(h4,'Color',[0.3,0.75 0.93])% Indigo
set(h12,'Color',[0.6,0.2,0])%brown
set(h2,'Color',[0.6,0.2,0])%brown
set(h3,'Color',[0.6,0.2,0])%brown

hl=legend([h1,h12,h2,h3,h4,h5],'DBF:WMMSE','PS-HBF:Ideal','PS-HBF:2bits','PS-HBF:1bits','SW-HBF:PGA-TS','SW-HBF:PGA-random');
set(hl,'FontSize',10,'FontWeight','bold')
ylabel('Average Weighted Sum Rate (bits/s/Hz)')
xlabel('Bandwidth (GHz)')
% xlabel('BSR')
xticks(xx)
xlim([xx(1) xx(end)])
grid on
box on

xx=B_all/GHz;
% xx=BSR_all;
figure
hold on
h1=plot(xx,EE_DBF,'k-*','LineWidth',1.5);
h12=plot(xx,EE_PS_4b,'--o','LineWidth',1.5);
h2=plot(xx,EE_PS_2b,'-o','LineWidth',1.5);
h3=plot(xx,EE_PS_1b,':o','LineWidth',1.5);
h4=plot(xx,EE_SW_pgats,'r-s','LineWidth',1.5);
h5=plot(xx,EE_SW_rdn,'b-^','LineWidth',1.5);

% set(h5,'Color',[0,0.5,0])%牛油绿
% set(h11,'Color',[0.4940 0.1840 0.5560])% purple
% set(h4,'Color',[0.3,0.75 0.93])% Indigo
set(h2,'Color',[0.6,0.2,0])%brown
set(h3,'Color',[0.6,0.2,0])%brown

hl=legend([h1,h12,h2,h3,h4,h5],'DBF:WMMSE','PS-HBF:Ideal','PS-HBF:2bits','PS-HBF:1bits','SW-HBF:PGA-TS','SW-HBF:PGA-random');
set(hl,'FontSize',10,'FontWeight','bold')
ylabel('Average EE (Gbits/J)')
xlabel('Bandwidth (GHz)')
% xlabel('BSR')
xticks(xx)
xlim([xx(1) xx(end)])
grid on
box on

%% Nu=8
xx=B_all/GHz;
% xx=BSR_all;
figure
hold on

h1=plot(xx,SE_DBF_avr_v1,'k-*','LineWidth',1.5);
h12=plot(xx,SE_PS_4b_avr_v1,'--o','LineWidth',1.5);
h2=plot(xx,SE_PS_2b_avr_v1,'-o','LineWidth',1.5);
h3=plot(xx,SE_PS_1b_avr_v1,':o','LineWidth',1.5);
h4=plot(xx,SW_PgaTS_avr_v1,'r-s','LineWidth',1.5);
h5=plot(xx,SW_PgaRdn_avr_v1,'b-^','LineWidth',1.5);
% set(h5,'Color',[0,0.5,0])%牛油绿
% set(h11,'Color',[0.4940 0.1840 0.5560])% purple
% set(h4,'Color',[0.3,0.75 0.93])% Indigo
set(h12,'Color',[0.6,0.2,0])%brown
set(h2,'Color',[0.6,0.2,0])%brown
set(h3,'Color',[0.6,0.2,0])%brown

hl=legend([h1,h12,h2,h3,h4,h5],'DBF:WMMSE','PS-HBF:Ideal','PS-HBF:2bits','PS-HBF:1bits','SW-HBF:PGA-TS','SW-HBF:PGA-random');
set(hl,'FontSize',10,'FontWeight','bold')
ylabel('Average Weighted Sum Rate (bits/s/Hz)')
xlabel('Bandwidth (GHz)')
% xlabel('BSR')
xticks(xx)
xlim([xx(1) xx(end)])
grid on
box on

xx=B_all/GHz;
% xx=BSR_all;
figure
hold on
h1=plot(xx,EE_DBF_v1,'k-*','LineWidth',1.5);
h12=plot(xx,EE_PS_4b_v1,'--o','LineWidth',1.5);
h2=plot(xx,EE_PS_2b_v1,'-o','LineWidth',1.5);
h3=plot(xx,EE_PS_1b_v1,':o','LineWidth',1.5);
h4=plot(xx,EE_SW_pgats_v1,'r-s','LineWidth',1.5);
h5=plot(xx,EE_SW_rdn_v1,'b-^','LineWidth',1.5);

% set(h5,'Color',[0,0.5,0])%牛油绿
% set(h11,'Color',[0.4940 0.1840 0.5560])% purple
% set(h4,'Color',[0.3,0.75 0.93])% Indigo
set(h2,'Color',[0.6,0.2,0])%brown
set(h3,'Color',[0.6,0.2,0])%brown

hl=legend([h1,h12,h2,h3,h4,h5],'DBF:WMMSE','PS-HBF:Ideal','PS-HBF:2bits','PS-HBF:1bits','SW-HBF:PGA-TS','SW-HBF:PGA-random');
set(hl,'FontSize',10,'FontWeight','bold')
ylabel('Average EE (Gbits/J)')
xlabel('Bandwidth (GHz)')
% xlabel('BSR')
xticks(xx)
xlim([xx(1) xx(end)])
grid on
box on

ccc=1;

end

