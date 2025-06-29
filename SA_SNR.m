function TC=SA_SNR
clc;clear all;close all
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

episilon=1e-4;

Pb_dBm=30;
Pb=10^(Pb_dBm/10);%mW single carrier

SNRx_dB=-10:5:30;
SNRx=10.^(SNRx_dB/10);
noise_power_all=Pb./SNRx;% noise power

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

tol=1e-4;T=1e+3;kappa=0.95;
params.T=T;
params.tol=tol;
params.kappa=kappa;

len=length(SNRx_dB);

BSR_score=AtnDis*Nt*Frac_Bw/8;
avr=100;
TC=zeros(1,avr);
%%
IniMatrix=zeros(avr,len);
SE_DBF=IniMatrix;

SW_rs=IniMatrix;
SW_Rdn_TS=IniMatrix; 
SW_Apga_TS=IniMatrix;
SW_SA=IniMatrix;
SW_random=IniMatrix; 


Flag_Identity_Ini=1;  
obj_fun_pg=@compute_rate_v0;        
obj_fun=@compute_rate_v0;
parfor n=1:avr
    [Hcc,ch_params]=Channel_Gen_compact(Nt,Nr,L,AtnDis,B,ofdm_paras);
    for id=1:len   
         noise_power=noise_power_all(id);
        [SE_DBF(n,id),~,DBF,Proj_Comb]=DBF_water_filling(Hcc,Pb,Ns,noise_power);
        
        SW_SA(n,id)=HBF_SA(params,Pb,Ns,N_rf,Hcc,noise_power,episilon); 
        
        [SW_Apga_TS(n,id),~]=HBF_APGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,noise_power,episilon,delta,obj_fun,obj_fun_pg,Flag_Identity_Ini,Proj_Comb);

        [SW_Rdn_TS(n,id),~]=HBF_TS(params_ts,Pb,Ns,N_rf,Hcc,noise_power,episilon,obj_fun,Flag_Identity_Ini,Proj_Comb); 

        [SW_rs(n,id),~]=HBF_APGA_TbRs(params_pga,params_TbRs,Pb,Ns,N_rf,Hcc,noise_power,episilon,delta,obj_fun,obj_fun_pg,Flag_Identity_Ini,Proj_Comb);

        [SW_random(n,id),~]=SW_HBF_random(Pb,Ns,N_rf,Hcc,noise_power,episilon);
    end
end
ccc=1
time_mark=datestr(now,'mmmm-dd');
file_name=['SA_SNR_',time_mark];
save(file_name,'-v7.3')
len_avr=avr;
% x_num=7;
SE_DBF_avr=mean(SE_DBF(1:len_avr,:));

SW_SA_avr=mean(SW_SA(1:len_avr,:)); 
SW_Apga_TS_avr=mean(SW_Apga_TS(1:len_avr,:)); 
SW_Rdn_TS_avr=mean(SW_Rdn_TS(1:len_avr,:)); 
SW_rs_avr=mean(SW_rs(1:len_avr,:));
SW_random_avr=mean(SW_random(1:len_avr,:)); 


x=SNRx_dB;
figure
hold on
h0=plot(x,SE_DBF_avr,'k-*','LineWidth',1.5);
h1=plot(x,SW_SA_avr,'b-s','LineWidth',1.5);
h2=plot(x,SW_Apga_TS_avr,'r-d','LineWidth',1.5);
h3=plot(x,SW_Rdn_TS_avr,'--v','LineWidth',1.5);
h4=plot(x,SW_rs_avr,':p','LineWidth',1.5);
h5=plot(x,SW_random_avr,'m-x','LineWidth',1.5);

set(h3,'Color',[0,0.5,0])%牛油绿
set(h4,'Color',[0.6,0.2,0])%brown

legend([h0,h1,h2,h3,h4,h5],'DBF','SA','PGA-TS','TS','PGA-TbRs','Random');
ylabel('Average SE (bits/s/Hz)')
xlabel('SNR(dB)')
% xticks(Nall)
xlim([SNRx_dB(1) SNRx_dB(end)])
grid on
box on

ccc=1;

end









