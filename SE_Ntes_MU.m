function TC=SE_Ntes_MU
clc;clear all;close all
tic
rng(0)
GHz=1e+9;%bandwidth=1GHz
Nc=64;% carrier number
D=Nc/4;% maximum delay
B=30*GHz;% bandwidth
fc=300*GHz;% carrier frequency
% Nt=64;
Nu=4;
% Nrf=Nu;
L=4;
AtnDis=1/2;
Frac_Bw=B/fc;

ofdm_paras.B=B;
ofdm_paras.Nc=Nc;% carrier number
ofdm_paras.D=D;% maximum delay
ofdm_paras.fc=fc;% carrier frequency

Pb_dBm=30;
Pb=10^(Pb_dBm/10);%mW single carrier

SNRx_dB=10;
% SNRx_dB=-10:5:30;
SNRx=10.^(SNRx_dB/10);
noise_power=Pb./SNRx;% noise power

%% Intialization
tol=1e-4;
params_pga.acc=1e-3;
params_pga.iter_max=1e+3;

Iter_max=200;
params_ts.iter_max=Iter_max;
params_ts.local_iter_max=10;%TS stop criteria
params_ts.Len_tl=Iter_max;% length of tabu list

% Wgtu=rand([1,Nu]);
Wgtu=ones(1,Nu);
iter_max=1e+3;

N_all=2.^(3:8);
len=length(N_all);
avr=100;

%% 
IniMatrix=zeros(avr,len);
DBF_all=IniMatrix;
ES_all1=IniMatrix;
PgaTs_all1=IniMatrix;
PgaRdn_all1=IniMatrix;  

ES_all2=IniMatrix;
PgaTs_all2=IniMatrix;
PgaRdn_all2=IniMatrix; 
parfor n=1:avr


    for id=1:len   
         Nt=N_all(id);
     
         Hall=Channel_MU(Nt,Nu,L,AtnDis,B,ofdm_paras);
%          DBF_all(n,id)=WMMSE_MISO(Hall,Wgtu,Pb,noise_power,50,tol);

         Nrf=4;
         [PgaTs_all1(n,id),ES_all1(n,id),PgaRdn_all1(n,id)]=TwoStepHBF_design_MISO(params_ts,params_pga,Hall,Wgtu,Nrf,Pb,noise_power,iter_max,tol);
        
%          Nrf=8;
%          [PgaTs_all2(n,id),ES_all2(n,id),PgaRdn_all2(n,id)]=TwoStepHBF_design_MISO(params_ts,params_pga,Hall,Wgtu,Nrf,Pb,noise_power,iter_max,tol);
                                   
        ccc=1;
    end
%     TC(n)=toc/60;
%     disp([' avr=',num2str(n)])

end
TC=toc;
time_mark=datestr(now,'mmmm-dd');
file_name=['SE_Ntes_MU_',time_mark];
save(file_name,'-v7.3')

len_avr=avr;
ES_avr1=mean(ES_all1(1:len_avr,:));
PgaTs1=mean(PgaTs_all1(1:len_avr,:));
PgaRdn1=mean(PgaRdn_all1(1:len_avr,:)); 

ES_avr2=mean(ES_all2(1:len_avr,:));
PgaTs2=mean(PgaTs_all2(1:len_avr,:));
PgaRdn2=mean(PgaRdn_all2(1:len_avr,:));


time_mark=datestr(now,'mmmm-dd');
file_name=['SE_Ntes_MU_',time_mark];
save(file_name,'-v7.3')
   
altstyles = {'m--o'; 'r-p';'g-h'; 'b-s'; 'c-v'; 'k-*'; 'r-.p';'g-.h'; 'b-.s'; 'c-.v'; 'k-.*'};
% x=B_all/GHz;
x=N_all;
figure
hold on
h1=plot(x,ES_avr1,'r-*',x,PgaTs1,'b-s',x,PgaRdn1,'k-o','LineWidth',1.5);
% h2=plot(x,ES_avr2,'r--*',x,PgaTs2,'b--s',x,PgaRdn2,'k--o','LineWidth',1.5);
legend([h1(1),h1(2),h1(3)],'ES','Algorithm 5','PGA-random');
% legend([h1(1),h1(2),h1(3),h2(1),h2(2),h2(3)],'ES $N_{\rm RF}=4$','Algorithm 5 $N_{\rm RF}=4$','PGA-random $N_{\rm RF}=4$', ...
%     'ES $N_{\rm RF}=8$','Algorithm 5 $N_{\rm RF}=8$','PGA-random $N_{\rm RF}=8$','interpreter','latex');
ylabel('Average Weighted Sum Rate (bits/s/Hz)')
xlabel('Number of transmit antennas ($N_{\rm T}$)','interpreter','latex')
xticks(x)
xlim([x(1) x(end)])
grid on
box on
ccc=1;