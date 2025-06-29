function TC=SE_Nt_es
clc;clear all;close all
%% variable: SNR
rng(0)
GHz=1e+9;%bandwidth=1GHz
Nc=64;% carrier number
D=Nc/4;% maximum delay
B=30*GHz;% bandwidth
fc=300*GHz;% carrier frequency
% Nt=16;
% Nr=Nt;
N_rf=2;
Ns=N_rf;
AtnDis=1/2;
Frac_Bw=B/fc;
% BSR_coef=AtnDis*Nt*Frac_Bw/8;
ofdm_paras.B=B;
ofdm_paras.Nc=Nc;% carrier number
ofdm_paras.D=D;% maximum delay
ofdm_paras.fc=fc;% carrier frequency
L=4;

episilon=1e-4;

Pb_dBm=30;
Pb=10^(Pb_dBm/10);%mW

SNRx_dB=10;
SNRx=10.^(SNRx_dB/10);
noise_power=Pb./SNRx;% noise power


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


N_all=4:2:12;
len=length(N_all);
avr=100;
TC=zeros(1,avr);
%% 
IniMatrix=zeros(avr,len);
SE_DBF=IniMatrix;

SE_PS_LSAA_4b=IniMatrix;
SE_PS_LSAA_4b_tx=IniMatrix;

SE_PS_LSAA_2b=IniMatrix;
SE_PS_LSAA_1b=IniMatrix;

SW_es=IniMatrix;
SW_rs=IniMatrix;
SW_random=IniMatrix;  
SW_Apga_TS=IniMatrix;  
SW_Rdn_TS=IniMatrix;  
 
Flag_Identity_Ini=1;  
obj_fun_pg=@compute_rate_v0;
obj_fun=@compute_rate_v0;  
parfor n=1:avr

% load('test_es.mat','Hcc');
%     tic
    for id=1:len   
         Nt=N_all(id);
         Nr=Nt;
         [Hcc,~]=Channel_Gen_compact(Nt,Nr,L,AtnDis,B, ofdm_paras)                
         [SE_DBF(n,id),~,DBF,Proj_Comb]=DBF_water_filling(Hcc,Pb,Ns,noise_power); 

         [SW_es(n,id),~]=SW_HBF_ES(Pb,Ns,N_rf,Hcc,noise_power,episilon,obj_fun);

         [SW_Apga_TS(n,id),~]=HBF_APGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,noise_power,episilon,delta,obj_fun,obj_fun_pg,Flag_Identity_Ini,Proj_Comb);
         [SW_Rdn_TS(n,id),~]=HBF_TS(params_ts,Pb,Ns,N_rf,Hcc,noise_power,episilon,obj_fun,Flag_Identity_Ini,Proj_Comb);      
         [SW_rs(n,id),~]=HBF_APGA_TbRs(params_pga,params_TbRs,Pb,Ns,N_rf,Hcc,noise_power,episilon,delta,obj_fun,obj_fun_pg,Flag_Identity_Ini,Proj_Comb);  
         [SW_random(n,id),~]=SW_HBF_random(Pb,Ns,N_rf,Hcc,noise_power,episilon);      
                               

        ccc=1;
    end
%     TC(n)=toc/60;
%     disp([' avr=',num2str(n)])

end
time_mark=datestr(now,'mmmm-dd');
file_name=['SE_Ntes_',time_mark];
save(file_name,'-v7.3')
len_avr=avr;
SW_es_avr=mean(SW_es(1:len_avr,:));
SW_rs_avr=mean(SW_rs(1:len_avr,:));
SW_random_avr=mean(SW_random(1:len_avr,:)); 
SW_Apga_TS_avr=mean(SW_Apga_TS(1:len_avr,:)); 
SW_Rdn_TS_avr=mean(SW_Rdn_TS(1:len_avr,:)); 


% time_mark=datestr(now,'mmmm-dd');
% file_name=['SE_SNR_es',time_mark];
%  save(file_name,'-v7.3')
   
altstyles = {'m--o'; 'r-p';'g-h'; 'b-s'; 'c-v'; 'k-*'; 'r-.p';'g-.h'; 'b-.s'; 'c-.v'; 'k-.*'};
% x=B_all/GHz;
x=N_all;
figure
hold on
h0=plot(x,SW_es_avr,'r-*','LineWidth',1.5);
h1=plot(x,SW_Apga_TS_avr,'b-d','LineWidth',1.5);
h2=plot(x,SW_Rdn_TS_avr,'k--s','LineWidth',1.5);
h3=plot(x,SW_rs_avr,':o','LineWidth',1.5);
h4=plot(x,SW_random_avr,'-x','LineWidth',1.5);

set(h3,'Color',[0,0.5,0])%牛油绿
set(h4,'Color',[0.6,0.2,0])%brown


legend([h0,h1,h2,h3,h4],'ES','SW-HBF:PGA-TS','SW-HBF:TS','SW-HBF:PGA-TbRs','SW-HBF:random');
ylabel('Average SE (bits/s/Hz)')
xlabel('Number of antennas')
% xticks(Nall)
xlim([x(1) x(end)])
grid on
box on


ccc=1;
end
