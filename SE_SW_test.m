clc;clear all;close all
G=1e+9;%bandwidth=1GHz
Nc=64;% carrier number
D=Nc/4;% maximum delay, see configure of simulation in 'Dynamic Subarrays for Hybrid Precoding inWideband mmWave MIMO Systems'
B=30*G;% bandwidth
fc=300*G;% carrier frequency
AtnDis=1/2;
Nt=64;
Nr=Nt;
N_rf=4;
Ns=N_rf;% number of data stream
Z=Nr/N_rf;
% Size_cob=factorial(Nr)/factorial(Z)^N_rf;
L=4;% path number

ofdm_paras.Nc=Nc;% carrier number
ofdm_paras.D=D;% maximum delay
ofdm_paras.fc=fc;% carrier frequency
ofdm_paras.B=B;

Pb_dBm=30;
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

[Hcc,~]=Channel_Gen_compact(Nt,Nr,L,AtnDis,B,ofdm_paras);%(Nr,Nt,Nc)

        % [SE_DBF,DBF]=DBF_water_filling(Hcc,Pb,Ns,noise_power);
        % 
        % 
        %  [SE_PS_pp_qt,W_rf_ps]=HBF_pp_quantized(Hcc,Pb,Ns,N_rf,noise_power,episilon,phase_Set);  
        %  SE_HBF=HBF_LSAA(Hcc,Pb,Ns,N_rf,noise_power,episilon,phase_Set);
        %  [SE_PS_pp2_qt,W_rf_ps2]=HBF_pp2_quantized(Hcc,Pb,Ns,N_rf,noise_power,episilon,phase_Set);  
        % for c=1:Nc
        %     Ini1.proj_comb{c}=eye(Nr);
        %     Ini2.proj_comb{c}=DBF.proj_comb{c};
        % end 
% 
%            
%  Nall=2.^(3:7);   
%  Ln=length( Nall);
%  MC=10;
%  Dall=zeros(MC,Ln);
%  Dpg=zeros(MC,Ln);
%  Dpgr=zeros(MC,Ln);
% for avr=1:MC
%     for id=1:Ln
%         Nt=Nall(id);
%         Nr=Nt;
%         for c=1:Nc
%             Ini1.proj_comb{c}=eye(Nr);
%             Ini2.proj_comb{c}=DBF.proj_comb{c};
%         end 
%            
%         [Hcc,~]=Channel_Gen_compact(Nt,Nr,L,AtnDis,ofdm_paras);%(Nr,Nt,Nc)
%         obj_fun_pg=@compute_rate_v0;
%         [Dall(avr,id),Dpg(avr,id), Dpgr(avr,id)]=dimension_reduction(params_pga,Hcc,Pb,Ns,N_rf,Ini1,noise_power,obj_fun_pg,delta);
%         disp(['id=',num2str(id),', avr=',num2str(avr)])
%     end
% end
% Dall_avr=mean(Dall);
% Dpg_avr=mean(Dpg);
% Dpgr_avr=mean( Dpgr);
% x= Nall;
% figure
% plot(x,Dall_avr,'r-o',x,1/2*Dall_avr,'r--o',x,Dpg_avr,'b-s',x,Dpgr_avr,'k-p')
% legend('Orginal','1/2 of Orginal','after PGA','after PGA+refinement')
% xlabel('Number of antennas')
% xticks(Nall)
% xlim([Nall(1) Nall(end)])
% ylabel('Dimension of feasible space')
% grid on
% box on


obj_fun_pg=@compute_rate_v0;
obj_fun=@compute_rate_v0;  
 Flag_Identity_Ini=1; 
[SE_DBF,~,DBF,Proj_Comb]=DBF_water_filling(Hcc,Pb,Ns,noise_power);

t1=tic;
NumNeighbor=1e+4;
[SW_Apga_TS1,W_rf_Apga,Tc_PgaTs]=HBF_APGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,noise_power,episilon,delta, ...
             obj_fun,obj_fun_pg,Flag_Identity_Ini,Proj_Comb,NumNeighbor);
tc1=toc(t1)

t2=tic;
NumNeighbor=8;
[SW_Apga_TS2,W_rf_Apga,Tc_PgaTs]=HBF_APGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,noise_power,episilon,delta, ...
             obj_fun,obj_fun_pg,Flag_Identity_Ini,Proj_Comb,NumNeighbor);
tc2=toc(t2)

t3=tic;
NumNeighbor=16;
[SW_Apga_TS3,W_rf_Apga,Tc_PgaTs]=HBF_APGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,noise_power,episilon,delta, ...
             obj_fun,obj_fun_pg,Flag_Identity_Ini,Proj_Comb,NumNeighbor);
tc3=toc(t3)
         % [SW_Rdn_TS,W_rf_TS, Tc_TS]=HBF_TS(params_ts,Pb,Ns,N_rf,Hcc,noise_power,episilon,obj_fun,Flag_Identity_Ini,Proj_Comb);  
         % t0_1=clock;
         % [SW_rdn_ts1_ini1,W_rf0_ini1,tc_TS]=HBF_TS(params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,obj_fun);
%          [SW_rdn_ts1_ini2,W_rf0_ini2]=HBF_TS(params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,obj_fun);

%          [SW_npag_ts_ini1,W_rf0_npg_ini1]=HBF_NPGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,obj_fun,obj_fun_pg);
%          [SW_npag_ts_ini2,W_rf0_npg_ini2]=HBF_NPGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,Ini2,noise_power,episilon,obj_fun,obj_fun_pg);
         
%         t1=clock;
%         obj_fun=@compute_rate_v1;
%         obj_fun_pg=@compute_rate_v0;
%          [SW_pga_ts2,W_rf1_pg]=HBF_APGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,delta,obj_fun,obj_fun_pg);
%          t1_1=clock;
%          [SW_rdn_ts2,W_rf1]=HBF_TS(params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,obj_fun);
% 
%          t2=clock;
% 
%          obj_fun=@compute_rate_v2;
%          obj_fun_pg=@compute_rate_v0;
%          [SW_pga_ts3,W_rf2_pg]=HBF_APGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,delta,obj_fun,obj_fun_pg);
%           t2_1=clock;
%          [SW_rdn_ts3,W_rf2]=HBF_TS(params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,obj_fun);
% 
%          t3=clock;
%           obj_fun=@compute_rate_v3;
%           obj_fun_pg=@compute_rate_v0;
%          [SW_pga_ts4,W_rf3_pg]=HBF_APGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,delta,obj_fun,obj_fun_pg);
%          t3_1=clock;
%          [SW_rdn_ts4,W_rf3]=HBF_TS(params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,obj_fun);
% 
%          t4=clock;

         % obj_fun=@compute_rate_v0;
         % obj_fun_pg=@compute_rate_v0;
         % [SW_rs,W_rf_rs]=HBF_APGA_TbRs(params_pga,params_TbRs,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,delta,obj_fun,obj_fun_pg);
         % 
         % t5=clock;
%          [SW_random,W_rf_rdn]=SW_HBF_random(Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,obj_fun);

load('PGA_data.mat','data_pga');
load('TS_data.mat','data_ts');
load('TbRs_data.mat','data_TbRs');
load('PGA_TS_data.mat','data_pag_ts');
iter_all=1:Iter_max;
x_pga=1:length(data_pga);
x_ts=1:length(data_ts);
x_TbRs=1:length(data_TbRs);
x_pga_ts=1:length(data_pag_ts);

x_max=max([x_pga(end),x_ts(100),x_TbRs(end),x_pga_ts(end)]);
figure
plot(x_pga,data_pga,'b-',x_ts,data_ts,'k-',x_pga_ts,data_pag_ts,'g-',x_TbRs,data_TbRs,'r-')
legend('PGA Algorithm','TS Algorithm', 'PGA-TS Algorithm' ,'PGA-TbRs Algorithm')
xlabel('Iterations')
ylabel(' SE (bits/s/Hz)')
xlim([1 x_max])
grid on
box on
 



cccc=1;