
clc;clear all;close all
tic
rng(0)
GHz=1e+9;%bandwidth=1GHz
Nc=32;% carrier number
D=Nc/4;% maximum delay
B=30*GHz;% bandwidth
fc=300*GHz;% carrier frequency
Nt=256;
Nu=4;
Nrf=4;

L=4;

AtnDis=1/2;
Frac_Bw=B/fc;

ofdm_paras.B=B;
ofdm_paras.Nc=Nc;% carrier number
ofdm_paras.D=D;% maximum delay
ofdm_paras.fc=fc;% carrier frequency

radius=50;
r_h=10;
large_fading=generate_Large_fading(Nu,radius,r_h);

Pb_dBm=30;
Pb=10^(Pb_dBm/10);%mW single carrier

SNRx_dB=20;
% SNRx_dB=-10:5:30;
SNRx=10.^(SNRx_dB/10);
noise_power=Pb./SNRx;% noise power

Nr=Nu;
% [H_all,ch_params]=Channel_Gen_compact(Nt,Nr,L,AtnDis,B,ofdm_paras);


Hall=channel_MISO(Nt,Nu,L,AtnDis,B,ofdm_paras,large_fading);%(Nt,Nu,Nc)
% Hall=channel_THz_MISO(Nt,Nu,AtnDis,B,ofdm_paras);
% Hall=Channel_MU(Nt,Nu,L,AtnDis,B,ofdm_paras);
% load('MU_channel.mat')

Phase_Bit=[1 2 4];
phase_bit_num=length(Phase_Bit);
Nq_set=2.^Phase_Bit;
for lb=1:phase_bit_num   
phase_Set{lb}=(-pi:2*pi/Nq_set(lb):pi);
end
PhaseSet_1bit=phase_Set{1};
PhaseSet_2bit=phase_Set{2};
PhaseSet_4bit=phase_Set{3};

%% Intialization

params_pga.acc=1e-3;
params_pga.iter_max=1e+3;

Iter_max=200;
params_ts.iter_max=Iter_max;
params_ts.local_iter_max=10;%TS stop criteria
params_ts.Len_tl=Iter_max;% length of tabu list

% Wgtu=rand([1,Nu]);
Wgtu=ones(1,Nu);
iter_max=1e+3;tol=1e-4;

tic
% [SE_pgats,SE_es,SE_rdn]=TwoStepHBF_design_MISO(params_ts,params_pga,Hall,Wgtu,Nrf,Pb,noise_power,iter_max,tol)
toc
% SE_rdn=SWHBF_MISO_random(params_pga,Hall,Wgtu,Nrf,Pb,noise_power,iter_max,tol)
% SE_es=SWHBF_MISO_ES(params_pga,Hall,Wgtu,Nrf,Pb,noise_power,iter_max,tol)
SE_opt=WMMSE_MISO(Hall,Wgtu,Pb,noise_power,50,tol)


% SE_PS_4b=PSHBF_MISO_WeiYu(Hall,Nrf,Pb,noise_power,PhaseSet_2bit,iter_max,tol)
SE_PS_2b=PSHBF_MISO_WeiYu(Hall,Wgtu,Nrf,Pb,noise_power,PhaseSet_2bit,iter_max,tol)
SE_PS_1b=PSHBF_MISO_WeiYu(Hall,Wgtu,Nrf,Pb,noise_power,PhaseSet_1bit,iter_max,tol)





ccc=1;
%% defined function




