function TC=SE_SNR_MU
clc;clear all;close all
%% variable: SNR
rng(0)
tic
GHz=1e+9;%bandwidth=1GHz
Nc=64;% carrier number
D=Nc/4;% maximum delay
B=30*GHz;% bandwidth
fc=300*GHz;% carrier frequency
Nt=256;
Nu=4;
Nrf=Nu;

AtnDis=1/2;
Frac_Bw=B/fc;

ofdm_paras.B=B;
ofdm_paras.Nc=Nc;% carrier number
ofdm_paras.D=D;% maximum delay
ofdm_paras.fc=fc;% carrier frequency
L=4;

radius=50;
r_h=10;
large_fading=generate_Large_fading(Nu,radius,r_h);

Pb_dBm=30;
Pb=10^(Pb_dBm/10);%mW

SNRx_dB=-10:5:20;
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

Wgtu=ones(1,Nu);iter_max=1e+3;tol=1e-4;

len=length(SNRx_dB);
avr=100;
BSR_score=AtnDis*Nt*Frac_Bw/8;
%% 
IniMatrix=zeros(avr,len);
SE_DBF=IniMatrix;

SE_PS_4b=IniMatrix;
SE_PS_2b=IniMatrix;
SE_PS_1b=IniMatrix;

SW_PgaTS=IniMatrix;  
SW_PgaRdn=IniMatrix;  

         
parfor n=1:avr
    Hall=Channel_MU(Nt,Nu,L,AtnDis,B,ofdm_paras);
    % Hall=channel_MISO(Nt,Nu,L,AtnDis,B,ofdm_paras,large_fading);

    for id=1:len   
        noise_power=noise_power_all(id);
        
        SE_DBF(n,id)=WMMSE_MISO(Hall,Wgtu,Pb,noise_power,50,tol);

        [SW_PgaTS(n,id),~,SW_PgaRdn(n,id)]=TwoStepHBF_design_MISO(params_ts,params_pga,Hall,Wgtu,Nrf,Pb,noise_power,iter_max,tol);

        % SE_PS_4b=PSHBF_MISO_WeiYu(Hall,Nrf,Pb,noise_power,PhaseSet_2bit,iter_max,tol)
        SE_PS_4b(n,id)=PSHBF_MISO_WeiYu(Hall,Wgtu,Nrf,Pb,noise_power,PhaseSet_4bit,iter_max,tol);
        SE_PS_2b(n,id)=PSHBF_MISO_WeiYu(Hall,Wgtu,Nrf,Pb,noise_power,PhaseSet_2bit,iter_max,tol);
        SE_PS_1b(n,id)=PSHBF_MISO_WeiYu(Hall,Wgtu,Nrf,Pb,noise_power,PhaseSet_1bit,iter_max,tol);        

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


time_mark=datestr(now,'mmmm-dd');
file_name=['SE_SNR_MU_',time_mark];
save(file_name,'-v7.3')
   
%% draw figure
xx=SNRx_dB;
figure
hold on
h1=plot(xx,SE_DBF_avr,'k-*','LineWidth',1.5);
h12=plot(xx,SE_PS_4b_avr,'--o','LineWidth',1.5);
h2=plot(xx,SE_PS_2b_avr,'-o','LineWidth',1.5);
h3=plot(xx,SE_PS_1b_avr,':o','LineWidth',1.5);
h4=plot(xx,SW_PgaTS_avr,'r-s','LineWidth',1.5);
h5=plot(xx,SW_PgaRdn_avr,'b-^','LineWidth',1.5);

% set(h5,'Color',[0,0.5,0])%牛油绿
% set(h11,'Color',[0.4940 0.1840 0.5560])% purple
% set(h4,'Color',[0.3,0.75 0.93])% Indigo
set(h12,'Color',[0.6,0.2,0])%brown
set(h2,'Color',[0.6,0.2,0])%brown
set(h3,'Color',[0.6,0.2,0])%brown

hl=legend([h1,h12,h2,h3,h4,h5],'DBF:WMMSE','PS-HBF:Ideal','PS-HBF:2bits','PS-HBF:1bits','SW-HBF:PGA-TS','SW-HBF:PGA-random');
set(hl,'FontSize',10,'FontWeight','bold')
ylabel('Average Weighted Sum Rate (bits/s/Hz)')
xlabel('SNR (dB)')
% xticks(Nall)
xlim([xx(1) xx(end)])
grid on
box on

ccc=1;

end

