clc;clear all;close all
%% variable: SNR
rng(0)
GHz=1e+9;%bandwidth=1GHz
Nc=64;% carrier number
D=Nc/4;% maximum delay
B=30*GHz;% bandwidth
fc=300*GHz;% carrier frequency
Nt=64;
Nr=64;
N_rf=4;
Ns=N_rf;
AtnDis=1/2;
Frac_Bw=B/fc;

ofdm_paras.B=B;
ofdm_paras.Nc=Nc;% carrier number
ofdm_paras.D=D;% maximum delay
ofdm_paras.fc=fc;% carrier frequency
L=4;


Pb_dBm=30;
Pb=10^(Pb_dBm/10);%mW

SNRx_dB=10;
SNRx=10.^(SNRx_dB/10);
noise_power=Pb./SNRx;% noise power

delta=0.1;
params_pga.acc=1e-4;
params_pga.iter_max=1e+3;%parameters for pga

Iter_max=200;
params_ts.iter_max=Iter_max;
params_ts.local_iter_max=10;%TS stop criteria
params_ts.Len_tl=Iter_max;% length of tabu list

params_TbRs.iter_max=Iter_max;
params_TbRs.local_iter_max=10;
params_TbRs.Len_tl=Iter_max;

BSR_score=AtnDis*Nt*Frac_Bw/8;

obj_fun_pg=@compute_rate_v0;
obj_fun=@compute_rate_v0;  
[H_all,ch_params]=Channel_Gen_compact(Nt,Nr,L,AtnDis,B,ofdm_paras);


%% PGA algorithm
Np=10;
power_avr=Pb/Ns;
F_ini=rand(Nt,N_rf);
F_ini_int = mapping(F_ini,Ns);
coef=power_avr/noise_power;
Ht=0;
for c=1:Nc
    Hc=H_all(:,:,c);
    HBF_tx.H_eff{c}=Hc'*eye(Nr)*Hc;
end
[SE_pg,F_pg,data_pga]=PGA_new(params_pga,F_ini_int,HBF_tx,coef,obj_fun_pg);

%% PGA-TS
if rank(F_pg,1e-4)<Ns
    F_pg=F_ini_int;
end
[F_rf_ini,Vec_pgr,Index_bf,Npart1]=TS_pga_params(F_pg,Ns,delta);
if length(Index_bf)>=1
      NumNeighbor=1e+4;
[SE_SW_tx,F_rf,data_pga_ts]=Tabu_search_pga(params_ts,F_rf_ini,Ns,Vec_pgr,Index_bf,Npart1,HBF_tx,coef,obj_fun,NumNeighbor);

    NumNeighbor=20;
[SE_SW_tx32,F_rf32,data_pga_ts32]=Tabu_search_pga(params_ts,F_rf_ini,Ns,Vec_pgr,Index_bf,Npart1,HBF_tx,coef,obj_fun,NumNeighbor);

    NumNeighbor=10;
[SE_SW_tx16,F_rf16,data_pga_ts16]=Tabu_search_pga(params_ts,F_rf_ini,Ns,Vec_pgr,Index_bf,Npart1,HBF_tx,coef,obj_fun,NumNeighbor);

    NumNeighbor=5;
[SE_SW_tx8,F_rf8,data_pga_ts8]=Tabu_search_pga(params_ts,F_rf_ini,Ns,Vec_pgr,Index_bf,Npart1,HBF_tx,coef,obj_fun,NumNeighbor);
% save('PGA_TS_data.mat','data_pag_ts');
end
figure
h0=plot(data_pga_ts,'LineWidth',1.5);hold on
h1=plot(data_pga_ts32,'LineWidth',1.5);hold on
h2=plot(data_pga_ts16,'LineWidth',1.5);hold on
h3=plot(data_pga_ts8,'LineWidth',1.5);hold on
legend([h0,h1,h2,h3],'$N_{\rm nb}=|\mathcal{S}_I|$','$N_{\rm nb}=24$',...
    '$N_{\rm nb}=16$','$N_{\rm nb}=8$','interpreter','latex')
xlabel('Iterations')
ylabel('SE [bits/s/Hz]')
grid on
box on
%% TS alg
Nu=4;
% Hfp=zeros(Nt,Nu,Nc);
% for nc=1:Nc
%     Hc=H_all(:,:,nc);
%     Hfp(:,:,nc)=Hc.';
% end
Hall=channel_MISO(Nt,Nu,L,AtnDis,B,ofdm_paras);%(Nt,Nu,Nc)
% [SE_tx,F_rf,data_ts]=Tabu_search(params_ts,F_ini_int,Ns,HBF_tx,coef,obj_fun);
Wgtu=ones(1,Nu);iter_max=1e+3;tol=1e-4;
[Frf,Fbb_all,data_FP_MU,data_pga_MU]=FP_relaxed_HBF(params_pga,Hall,Wgtu,N_rf,Pb,noise_power,iter_max,tol);
%% PGA-TbRs alg

% Search_flag=pga_solution(F_pg,Ns);
% if Search_flag==1
%     F_pgr = rounding(F_pg,delta);
%     if rank(F_pgr,1e-4)<Ns
%        F_pgr= F_pg;
%     end
%     [vBest,F_rf,data_TbRs]=Tabu_RS(params_TbRs,F_pgr,Ns,HBF_tx,coef,obj_fun);
% %     save('TbRs_data.mat','data_TbRs');
% 
% end

time_mark=datestr(now,'mmmm-dd');
file_name=['Converge_data_',time_mark];
save(file_name,'-v7.3')
%% draw figure
data_pga=data_pga(1:80);
iter_all=1:Iter_max;
x_pga=1:length(data_pga);
x_fp=1:length(data_FP_MU);
% x_TbRs=1:length(data_TbRs);
x_pga_ts=1:length(data_pga_ts);

x_max=max([x_pga(end),x_fp(end),x_pga_ts(end)]);
figure
hold on
% plot(x_pga,data_pga,'b-',x_fp,data_ts,'k-',x_pga_ts,data_pag_ts,'g-',x_TbRs,data_TbRs,'r-','LineWidth',1.5)
% legend('PGA Algorithm','TS Algorithm', 'PGA-TS Algorithm' ,'PGA-TbRs Algorithm')
h1=plot(x_pga,data_pga/max(data_pga),'b-','LineWidth',1.5)
h2=plot(x_pga_ts,data_pga_ts/max(data_pga_ts),'r-','LineWidth',1.5)
h3=plot(x_fp,data_FP_MU/max(data_FP_MU),'k-','LineWidth',1.5)

legend([h1,h2,h3],'PGA Algorithm (SU-MIMO)','PGA-TS Algorithm (SU-MIMO)','FP-based Algorihtm (MU-MISO)' )
xlabel('Iterations')
ylabel('Normalized objective')
xlim([1 x_max])
grid on
box on
ccc=1;
%% defined function
function [W_rf_ini,Vec_pgr,Index_bf,Npart1]=TS_pga_params(W_pg,Ns,delta)
[N,N_rf]=size(W_pg);
W_pgr = rounding(W_pg,delta);
rank_pgr=rank(W_pgr,1e-4);
if rank_pgr<Ns
        W_pgr=W_pg;        
end
    
% W_pgr
% W_pg
Vec_pgr=vec(W_pgr);
W_rf_ini = mapping(W_pgr,Ns);
    
V_rf_ini=vec(W_rf_ini);
Diff_bf=abs(Vec_pgr-V_rf_ini);
[Index_bf,~]=find(vec(Diff_bf)>0);% Index_bf denotes the index in Vec_pgr

Determined_ele=setdiff((1:N*N_rf),Index_bf);
Vec_pgr_determined=Vec_pgr(Determined_ele);
Npart1=sum(Vec_pgr_determined);
end

function Search_flag=pga_solution(W_pg,Ns)
Vec_pg=vec(W_pg);
W_rf_ini = mapping(W_pg,Ns);
V_rf_ini=vec(W_rf_ini);
Diff_bf=abs(Vec_pg-V_rf_ini);
[Index_bf,~]=find(vec(Diff_bf)>0);% Index_bf denotes the index in Vec_pgr

    if length(Index_bf)>=1
        Search_flag=1;

    else
        Search_flag=0;
    end

end
