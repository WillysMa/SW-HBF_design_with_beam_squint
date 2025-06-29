function TC=SE_Nu_MU
clc;clear all;close all
%% variable: Number of User
rng(0)
tic
GHz=1e+9;%bandwidth=1GHz
Nc=64;% carrier number
D=Nc/4;% maximum delay
B=30*GHz;% bandwidth
fc=300*GHz;% carrier frequency
Nt=256;
% Nu=8;

AtnDis=1/2;

BSR_coef=AtnDis*Nt/fc/8;
ofdm_paras.B=B;
ofdm_paras.Nc=Nc;% carrier number
ofdm_paras.D=D;% maximum delay
ofdm_paras.fc=fc;% carrier frequency
L=4;

radius=50;
r_h=10;
% large_fading=generate_Large_fading(Nu,radius,r_h);

Pb_dBm=30;
Pb=10^(Pb_dBm/10);%mW

SNRx_dB=[-10 5 20];
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

params_pga.acc=1e-4;
params_pga.iter_max=1e+3;%parameters for pga

Iter_max=200;
params_ts.iter_max=Iter_max;
params_ts.local_iter_max=10;%TS stop criteria
params_ts.Len_tl=Iter_max;% length of tabu list

iter_max=1e+3;tol=1e-4;

Nu_all=1:8;
len=length(Nu_all);

avr=100;

%% 
IniMatrix=zeros(avr,len);
SE_DBF=IniMatrix;
SE_PS_4b=IniMatrix;
SE_PS_2b=IniMatrix;
SE_PS_1b=IniMatrix;
SW_PgaTS=IniMatrix;  
SW_PgaRdn=IniMatrix;  

SE_DBFl=IniMatrix;
SE_PS_4bl=IniMatrix;
SE_PS_2bl=IniMatrix;
SE_PS_1bl=IniMatrix;
SW_PgaTSl=IniMatrix;  
SW_PgaRdnl=IniMatrix;  

SE_DBFh=IniMatrix;
SE_PS_4bh=IniMatrix;
SE_PS_2bh=IniMatrix;
SE_PS_1bh=IniMatrix;
SW_PgaTSh=IniMatrix;  
SW_PgaRdnh=IniMatrix;  
         
parfor n=1:avr   
    for id=1:len   
        Nu=Nu_all(id);
        Nrf=Nu;
        Wgtu=ones(1,Nu);
        Hall=Channel_MU(Nt,Nu,L,AtnDis,B,ofdm_paras);
        % Hall=channel_MISO(Nt,Nu,L,AtnDis,B,ofdm_paras,large_fading);
        
        noise_power=noise_power_all(1);
        SE_DBFl(n,id)=WMMSE_MISO(Hall,Wgtu,Pb,noise_power,50,tol);
        [SW_PgaTSl(n,id),~,SW_PgaRdnl(n,id)]=TwoStepHBF_design_MISO(params_ts,params_pga,Hall,Wgtu,Nrf,Pb,noise_power,iter_max,tol);
        % SE_PS_4b=PSHBF_MISO_WeiYu(Hall,Nrf,Pb,noise_power,PhaseSet_2bit,iter_max,tol)
        SE_PS_4bl(n,id)=PSHBF_MISO_WeiYu(Hall,Wgtu,Nrf,Pb,noise_power,PhaseSet_4bit,iter_max,tol);
        SE_PS_2bl(n,id)=PSHBF_MISO_WeiYu(Hall,Wgtu,Nrf,Pb,noise_power,PhaseSet_2bit,iter_max,tol);
        SE_PS_1bl(n,id)=PSHBF_MISO_WeiYu(Hall,Wgtu,Nrf,Pb,noise_power,PhaseSet_1bit,iter_max,tol);        

        noise_power=noise_power_all(2);
        SE_DBF(n,id)=WMMSE_MISO(Hall,Wgtu,Pb,noise_power,50,tol);
        [SW_PgaTS(n,id),~,SW_PgaRdn(n,id)]=TwoStepHBF_design_MISO(params_ts,params_pga,Hall,Wgtu,Nrf,Pb,noise_power,iter_max,tol);
        % SE_PS_4b=PSHBF_MISO_WeiYu(Hall,Nrf,Pb,noise_power,PhaseSet_2bit,iter_max,tol)
        SE_PS_4b(n,id)=PSHBF_MISO_WeiYu(Hall,Wgtu,Nrf,Pb,noise_power,PhaseSet_4bit,iter_max,tol);
        SE_PS_2b(n,id)=PSHBF_MISO_WeiYu(Hall,Wgtu,Nrf,Pb,noise_power,PhaseSet_2bit,iter_max,tol);
        SE_PS_1b(n,id)=PSHBF_MISO_WeiYu(Hall,Wgtu,Nrf,Pb,noise_power,PhaseSet_1bit,iter_max,tol);  

        noise_power=noise_power_all(3);
        SE_DBFh(n,id)=WMMSE_MISO(Hall,Wgtu,Pb,noise_power,50,tol);
        [SW_PgaTSh(n,id),~,SW_PgaRdnh(n,id)]=TwoStepHBF_design_MISO(params_ts,params_pga,Hall,Wgtu,Nrf,Pb,noise_power,iter_max,tol);
        % SE_PS_4b=PSHBF_MISO_WeiYu(Hall,Nrf,Pb,noise_power,PhaseSet_2bit,iter_max,tol)
        SE_PS_4bh(n,id)=PSHBF_MISO_WeiYu(Hall,Wgtu,Nrf,Pb,noise_power,PhaseSet_4bit,iter_max,tol);
        SE_PS_2bh(n,id)=PSHBF_MISO_WeiYu(Hall,Wgtu,Nrf,Pb,noise_power,PhaseSet_2bit,iter_max,tol);
        SE_PS_1bh(n,id)=PSHBF_MISO_WeiYu(Hall,Wgtu,Nrf,Pb,noise_power,PhaseSet_1bit,iter_max,tol); 
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
SE_DBF_avr=mean(SE_DBF(1:len_avr,:));
SE_PS_1b_avr=mean(SE_PS_1b(1:len_avr,:));
SE_PS_2b_avr=mean(SE_PS_2b(1:len_avr,:));
SE_PS_4b_avr=mean(SE_PS_4b(1:len_avr,:));
SW_PgaTS_avr=mean(SW_PgaTS(1:len_avr,:));
SW_PgaRdn_avr=mean(SW_PgaRdn(1:len_avr,:)); 

SE_DBFl_avr=mean(SE_DBFl(1:len_avr,:));
SE_PS_1bl_avr=mean(SE_PS_1bl(1:len_avr,:));
SE_PS_2bl_avr=mean(SE_PS_2bl(1:len_avr,:));
SE_PS_4bl_avr=mean(SE_PS_4bl(1:len_avr,:));
SW_PgaTSl_avr=mean(SW_PgaTSl(1:len_avr,:));
SW_PgaRdnl_avr=mean(SW_PgaRdnl(1:len_avr,:)); 

SE_DBFh_avr=mean(SE_DBFh(1:len_avr,:));
SE_PS_1bh_avr=mean(SE_PS_1bh(1:len_avr,:));
SE_PS_2bh_avr=mean(SE_PS_2bh(1:len_avr,:));
SE_PS_4bh_avr=mean(SE_PS_4bh(1:len_avr,:));
SW_PgaTSh_avr=mean(SW_PgaTSh(1:len_avr,:));
SW_PgaRdnh_avr=mean(SW_PgaRdnh(1:len_avr,:)); 


time_mark=datestr(now,'mmmm-dd');
file_name=['SE_Nu_MU_',time_mark];
save(file_name,'-v7.3')
   
%% draw figure
lineSpec1={'r-o';'k-s';'b-+';'m-^';'g-d';'c-x';{[0.6,0.2,0],'-*'}};
lineSpec2={'r--o';'k--s';'b--+';'m--^';'g--d';'c--x';{[0.6,0.2,0],'--*'}};
lineSpec3={'r:o';'k:s';'b:+';'m:^';'g:d';'c:x';{[0.6,0.2,0],':*'}};
% lineSpec4={'r-.o';'k-.s';'b-.+';'m-.^';'g-.d';'c-.x';{[0.6,0.2,0],'-.*'}};
lineSpec0={'ro';'ks';'b+';'m^';'gd';'cx';{[0.6,0.2,0],'*'}};

lenxx=8;
xx=Nu_all(1:lenxx);
figure
hold on

lineSpec=lineSpec1;
% h11=plot(xx,SE_DBFl_avr(1:lenxx),lineSpec{1},'LineWidth',1.5);
h12=plot(xx,SE_PS_4bl_avr(1:lenxx),lineSpec{2},'LineWidth',1.5);
h22=plot(xx,SE_PS_2bl_avr(1:lenxx),lineSpec{3},'LineWidth',1.5);
h33=plot(xx,SE_PS_1bl_avr(1:lenxx),lineSpec{4},'LineWidth',1.5);
h44=plot(xx,SW_PgaTSl_avr(1:lenxx),lineSpec{1},'LineWidth',1.5);
h55=plot(xx,SW_PgaRdnl_avr(1:lenxx),lineSpec{5},'LineWidth',1.5);


lineSpec=lineSpec2;
% h11=plot(xx,SE_DBF_avr(1:lenxx),lineSpec{1},'LineWidth',1.5);
h12=plot(xx,SE_PS_4b_avr(1:lenxx),lineSpec{2},'LineWidth',1.5);
h22=plot(xx,SE_PS_2b_avr(1:lenxx),lineSpec{3},'LineWidth',1.5);
h33=plot(xx,SE_PS_1b_avr(1:lenxx),lineSpec{4},'LineWidth',1.5);
h44=plot(xx,SW_PgaTS_avr(1:lenxx),lineSpec{1},'LineWidth',1.5);
h55=plot(xx,SW_PgaRdn_avr(1:lenxx),lineSpec{5},'LineWidth',1.5);

lineSpec=lineSpec3;
% h11=plot(xx,SE_DBFh_avr(1:lenxx),lineSpec{1},'LineWidth',1.5);
h12=plot(xx,SE_PS_4bh_avr(1:lenxx),lineSpec{2},'LineWidth',1.5);
h22=plot(xx,SE_PS_2bh_avr(1:lenxx),lineSpec{3},'LineWidth',1.5);
h33=plot(xx,SE_PS_1bh_avr(1:lenxx),lineSpec{4},'LineWidth',1.5);
h44=plot(xx,SW_PgaTSh_avr(1:lenxx),lineSpec{1},'LineWidth',1.5);
h55=plot(xx,SW_PgaRdnh_avr(1:lenxx),lineSpec{5},'LineWidth',1.5);

lineSpec=lineSpec0;
% h1=plot(xx,SE_DBF_avr(1:lenxx),lineSpec{1},'LineWidth',1.5);
h12=plot(xx,SE_PS_4b_avr(1:lenxx),lineSpec{2},'LineWidth',1.5);
h2=plot(xx,SE_PS_2b_avr(1:lenxx),lineSpec{3},'LineWidth',1.5);
h3=plot(xx,SE_PS_1b_avr(1:lenxx),lineSpec{4},'LineWidth',1.5);
h4=plot(xx,SW_PgaTS_avr(1:lenxx),lineSpec{1},'LineWidth',1.5);
h5=plot(xx,SW_PgaRdn_avr(1:lenxx),lineSpec{5},'LineWidth',1.5);
% set(h5,'Color',[0,0.5,0])%牛油绿
% set(h11,'Color',[0.4940 0.1840 0.5560])% purple
% set(h4,'Color',[0.3,0.75 0.93])% Indigo
% set(h2,'Color',[0.6,0.2,0])%brown
% set(h3,'Color',[0.6,0.2,0])%brown

% hl=legend([h1,h12,h2,h3,h4,h5],'DBF:WMMSE','PS-HBF:Ideal','PS-HBF:2bits','PS-HBF:1bits','SW-HBF:PGA-TS','SW-HBF:PGA-random');
hl=legend([h12,h2,h3,h4,h5],'PS-HBF:Ideal','PS-HBF:2bits','PS-HBF:1bits','SW-HBF:PGA-TS','SW-HBF:PGA-random');
set(hl,'FontSize',10,'FontWeight','bold')

ylabel('Average Weighted Sum Rate (bits/s/Hz)')
xlabel('Number of users ($N_{\rm U}$)','Interpreter','latex')
xticks(xx)
xlim([xx(1) xx(end)])
grid on
box on

ccc=1;

end

