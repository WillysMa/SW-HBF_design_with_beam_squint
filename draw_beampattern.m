clc;clear all;close all
G=1e+9;%bandwidth=1GHz
Nc=64;% carrier number
D=Nc/4;% maximum delay, see configure of simulation in 'Dynamic Subarrays for Hybrid Precoding inWideband mmWave MIMO Systems'
B=30*G;% bandwidth
fc=300*G;% carrier frequency
AtnDis=1/2;
Nt=64;
Nr=Nt;
N_rf=1;
Ns=N_rf;% number of data stream
Z=Nr/N_rf;
% Size_cob=factorial(Nr)/factorial(Z)^N_rf;
L=4;% path number

ofdm_paras.Nc=Nc;% carrier number
ofdm_paras.D=D;% maximum delay
ofdm_paras.fc=fc;% carrier frequency
ofdm_paras.B=B;

Pb_dBm=20;
Pb=10^(Pb_dBm/10);%mW
SNRx_dB=10;
SNRx=10.^(SNRx_dB/10);
noise_power=Pb./SNRx;% noise power
len=length(SNRx_dB);

Phase_Bit=4;
phase_bit_num=length(Phase_Bit);
Nq_set=2.^Phase_Bit;
phase_Set=(-pi:2*pi/Nq_set:pi);

% rng('default')
episilon=1e-6;
delta=0.1;
params_pga.acc=1e-4;
params_pga.iter_max=300;%parameters for pga

Iter_max=300;
params_ts.iter_max=Iter_max;
params_ts.local_iter_max=10;%TS stop criteria
params_ts.Len_tl=Iter_max;% length of tabu list

params_TbRs.iter_max=Iter_max;
params_TbRs.local_iter_max=10;
params_TbRs.Len_tl=Iter_max;

F_ini=rand(Nt,N_rf);
F_rf = mapping(F_ini,Ns);
Pf=F_rf*pinv(F_rf);
Vpf=eig(Pf);
Flag_Identity_Ini=1;
[Hcc,~]=Channel_Gen_compact(Nt,Nr,L,AtnDis,B,ofdm_paras);%(Nr,Nt,Nc)

[SE_DBF,DBF,~,Proj_Comb]=DBF_water_filling(Hcc,Pb,Ns,noise_power);


% [SE_PS_pp_qt,W_rf_ps]=HBF_pp_quantized(Hcc,Pb,Ns,N_rf,noise_power,episilon,phase_Set);
% % [SE_PS_LSAA_4b(n,id),~]=HBF_LSAA(Hcc,Pb,Ns,N_rf,noise_power,episilon,PhaseSet_4bit); 
% 
% receive_bf=W_rf_ps;
% beampattern(receive_bf,Ns,AtnDis,ofdm_paras)
% % SE_HBF=HBF_LSAA(Hcc,Pb,Ns,N_rf,noise_power,episilon,phase_Set);
% % [SE_PS_pp2_qt,W_rf_ps2]=HBF_pp2_quantized(Hcc,Pb,Ns,N_rf,noise_power,episilon,phase_Set);  

obj_fun=@compute_rate_v0;
obj_fun_pg=@compute_rate_v0;
[SW_pga_ts1_ini1,W_rf0_pg_ini1]=HBF_APGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,noise_power,episilon,delta,obj_fun,obj_fun_pg,Flag_Identity_Ini,Proj_Comb);

receive_bf=W_rf0_pg_ini1;
sum(receive_bf)
% beampattern(receive_bf,Ns,AtnDis,ofdm_paras)


% Define the variables
K = Nc; % Number of iterations for the outer summation
N = Nr;  % Number of iterations for the inner summation
b = B/fc;   % Value of b
Delta = AtnDis; % Value of Delta

% Generate binary values for w_n (0 or 1)
% w_n = randi([0, 1], 1, N); % Example: random binary vector of length N
w_n = receive_bf; % Example: random binary vector of length N
final_result=SumAbsSumCpxExp(w_n,K,N,Delta,b)

w_n=ones(1,Nr);
final_result=SumAbsSumCpxExp(w_n,K,N,Delta,b)
ccc=1;

function final_result=SumAbsSumCpxExp(w_n,K,N,Delta,b)

% Calculate the norm of w
norm_w = norm(w_n, 1);

% Initialize summation variable
sum_result = 0;

% Perform outer summation
for k = 1:K
    % Initialize inner summation variable
    inner_sum = 0;
    
    % Perform inner summation
    for n = 1:N
        % Compute complex exponential term
        exp_term = exp(-1j * 2 * Delta * pi * (1 + ((k / K) - ((K + 1) / (2 * K))) * b) * (n - 1));
        
        % Compute inner summation
        inner_sum = inner_sum + w_n(n) * exp_term;
    end
    
    % Compute absolute value of inner summation and add to outer summation
    sum_result = sum_result + abs(inner_sum);
end

% Compute final result
final_result = sum_result / (K * sqrt(N * norm_w));

end


