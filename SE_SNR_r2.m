function TC=SE_SNR_r2
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


Pb_dBm=30;
Pb=10^(Pb_dBm/10);%mW single carrier

SNRx_dB=-10:5:30;
SNRx=10.^(SNRx_dB/10);
noise_power_all=Pb./SNRx;% noise power

episilon=1e-6;
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
avr=50;
TC=zeros(1,avr);
%% 
IniMatrix=zeros(avr,len);
SE_DBF=IniMatrix;

SW_Apga_TS1=IniMatrix;
SW_Apga_TS2=IniMatrix;
SW_Apga_TS3=IniMatrix;


obj_fun_pg=@compute_rate_v0;
obj_fun=@compute_rate_v0;     
Flag_Identity_Ini=1;
tic
parfor n=1:avr
[Hcc,ch_params]=Channel_Gen_compact(Nt,Nr,L,AtnDis,B,ofdm_paras);

%     tic
    for id=1:len   
        noise_power=noise_power_all(id);
        [SE_DBF(n,id),~,DBF,Proj_Comb]=DBF_water_filling(Hcc,Pb,Ns,noise_power);


        NumNeighbor=1e+4;
        [SW_Apga_TS1(n,id),W_rf_Apga,Tc_PgaTs]=HBF_APGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,noise_power,episilon,delta, ...
                     obj_fun,obj_fun_pg,Flag_Identity_Ini,Proj_Comb,NumNeighbor);
             
        NumNeighbor=8;
        [SW_Apga_TS2(n,id),W_rf_Apga,Tc_PgaTs]=HBF_APGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,noise_power,episilon,delta, ...
                     obj_fun,obj_fun_pg,Flag_Identity_Ini,Proj_Comb,NumNeighbor);
        
        NumNeighbor=16;
        [SW_Apga_TS3(n,id),W_rf_Apga,Tc_PgaTs]=HBF_APGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,noise_power,episilon,delta, ...
                     obj_fun,obj_fun_pg,Flag_Identity_Ini,Proj_Comb,NumNeighbor);


        ccc=1;
    end
%     TC(n)=toc/60;
%     disp(['avr=',num2str(n)])
end
TC=toc

len_avr=avr;
% x_num=7;
SE_DBF_avr=mean(SE_DBF(1:len_avr,:));

SW_Apga_TS1_avr=mean(SW_Apga_TS1(1:len_avr,:)); 
SW_Apga_TS2_avr=mean(SW_Apga_TS2(1:len_avr,:)); 
SW_Apga_TS3_avr=mean(SW_Apga_TS3(1:len_avr,:)); 

time_mark=datestr(now,'mmmm-dd');
file_name=['SE_SNRr2','_',time_mark];
save(file_name,'-v7.3')
%% 
   
altstyles = {'m--o'; 'r-p';'g-h'; 'b-s'; 'c-v'; 'k-*'; 'r-.p';'g-.h'; 'b-.s'; 'c-.v'; 'k-.*'};
% x=B_all/GHz;
x=SNRx_dB;
figure
hold on
h0=plot(x,SW_Apga_TS1_avr,'r-*','LineWidth',1.5);
h1=plot(x,SW_Apga_TS2_avr,'b-v','LineWidth',1.5);
h2=plot(x,SW_Apga_TS3_avr,'k-o','LineWidth',1.5);


% set(h11,'Color',[0.4940 0.1840 0.5560])% purple

legend([h0,h2,h1],'$n_{\rm neighbor}=|\mathcal{S}_I|$','$n_{\rm neighbor}=16$','$n_{\rm neighbor}=8$','interpreter','latex');
ylabel('Average SE [bits/s/Hz]')
xlabel('SNR [dB]')
% xticks(Nall)
xlim([SNRx_dB(1) SNRx_dB(end)])
grid on
box on

% axesNew = axes('position',get(gca,'position'),'visible','off');
% legend(axesNew,[h8,h9,h10,h11],'SW-HBF:PGA-TS','SW-HBF:TS','SW-HBF:PGA-TbRs','SW-HBF:random','location','west')




ccc=1;


end



