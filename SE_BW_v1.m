function [SE_DBF_avr,SE_FTTD_avr]=SE_BW_v1(Ndt)
clc;clear all;close all
%% variable: Bandwidth
GHz=1e+9;%bandwidth=1GHz
Nc=64;% carrier number
D=Nc/4;% maximum delay
% B=10*GHz;% bandwidth
fc=300*GHz;% carrier frequency
Nt=512;
Nr=Nt;
N_rf=4;
Ns=N_rf;
AtnDis=1/2;

BSR_coef=AtnDis*Nt/fc/8;
% ofdm_paras.B=B;
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


Phase_Bit=[1 2 4];
phase_bit_num=length(Phase_Bit);
Nq_set=2.^Phase_Bit;
for lb=1:phase_bit_num   
phase_Set{lb}=(-pi:2*pi/Nq_set(lb):pi);
end
PhaseSet_1bit=phase_Set{1};
PhaseSet_2bit=phase_Set{2};
PhaseSet_4bit=phase_Set{3};

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



BSR_all=0.1:0.5:3.5;
B_all=BSR_all/BSR_coef;
len=length(B_all);

% B_all/GHz

N_ttd0=ceil(Nt/2*max(B_all)/fc);
N_ttd1=N_ttd0;
while mod(Nt,N_ttd1)
    N_ttd1=N_ttd1+1;
end
% Ndt=2;
Ndr=Ndt;
Nfdt=Nt/N_rf;
Nfdr=Nr/N_rf;

saved_name='SE_BWttd'+num2str(Ndt);
avr=10;

%% 
IniMatrix=zeros(avr,len);
SE_DBF=IniMatrix;
SE_FTTD=IniMatrix;
SE_TTD_4b=IniMatrix;
SE_TTD_2b=IniMatrix;
SE_TTD_1b=IniMatrix;

SE_PS_LSAA_4b=IniMatrix;
SE_PS_LSAA_2b=IniMatrix;
SE_PS_LSAA_1b=IniMatrix;

SW_rs=IniMatrix;
SW_random=IniMatrix;  
SW_Apga_TS=IniMatrix;  
SW_Rdn_TS=IniMatrix;  

for c=1:Nc
    Ini1.proj_comb{c}=eye(Nr);
end 
obj_fun_pg=@compute_rate_v0;        

for n=1:avr

    parfor id=1:len 
        B=B_all(id);
%         ofdm_paras.B=B_all(id);
        [Hcc,ch_params]=Channel_Gen_compact(Nt,Nr,L,AtnDis,B,ofdm_paras);
        
        [SE_DBF(n,id),~,DBF]=DBF_water_filling(Hcc,Pb,Ns,noise_power);
%          [SE_PS_pp_1b(n,id),~]=HBF_pp_quantized(Hcc,Pb,Ns,N_rf,noise_power,episilon,phase_Set{1}); 
%          [SE_PS_pp_2b(n,id),~]=HBF_pp_quantized(Hcc,Pb,Ns,N_rf,noise_power,episilon,phase_Set{2});   
%          [SE_PS_pp_4b(n,id),~]=HBF_pp_quantized(Hcc,Pb,Ns,N_rf,noise_power,episilon,phase_Set{3});  

         [SE_PS_LSAA_1b(n,id),~]=HBF_LSAA(Hcc,Pb,Ns,N_rf,noise_power,episilon,PhaseSet_1bit); 
         [SE_PS_LSAA_2b(n,id),~]=HBF_LSAA(Hcc,Pb,Ns,N_rf,noise_power,episilon,PhaseSet_2bit);   
         [SE_PS_LSAA_4b(n,id),~]=HBF_LSAA(Hcc,Pb,Ns,N_rf,noise_power,episilon,PhaseSet_4bit); 

         [SE_TTD_1b(n,id),~]=HBF_DDP(Hcc,Pb,Ns,N_rf,Ndt,Ndr,noise_power,AtnDis,B,ch_params,ofdm_paras,PhaseSet_1bit);
         [SE_TTD_2b(n,id),~]=HBF_DDP(Hcc,Pb,Ns,N_rf,Ndt,Ndr,noise_power,AtnDis,B,ch_params,ofdm_paras,PhaseSet_2bit);
         [SE_TTD_4b(n,id),~]=HBF_DDP(Hcc,Pb,Ns,N_rf,Ndt,Ndr,noise_power,AtnDis,B,ch_params,ofdm_paras,PhaseSet_4bit);

         [SE_FTTD(n,id),~]=HBF_FTTD(Hcc,DBF,N_rf,Nfdt,Nfdr,noise_power,AtnDis,B,ofdm_paras);
         
         obj_fun=@compute_rate_v0;        
         [SW_Apga_TS(n,id),~]=HBF_APGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,delta,obj_fun,obj_fun_pg);
         [SW_Rdn_TS(n,id),~]=HBF_TS(params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,obj_fun);      

         [SW_rs(n,id),~]=HBF_APGA_TbRs(params_pga,params_TbRs,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,delta,obj_fun,obj_fun_pg);  
         [SW_random(n,id),~]=SW_HBF_random(Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,obj_fun);





        ccc=1;
    end
    disp([' avr=',num2str(n)])
    time_mark=datestr(now,'mmmm-dd');
    file_name=[saved_name,'_',time_mark];
    save(file_name,'-v7.3')
end

len_avr=12;
x_num=8;

SE_DBF_avr=mean(SE_DBF(1:len_avr,:));
SE_FTTD_avr=mean(SE_FTTD(1:len_avr,:));
SE_TTD_1b_avr=mean(SE_TTD_1b(1:len_avr,:));
SE_TTD_2b_avr=mean(SE_TTD_2b(1:len_avr,:));
SE_TTD_4b_avr=mean(SE_TTD_4b(1:len_avr,:));

SE_PS_LSAA_1b_avr=mean(SE_PS_LSAA_1b(1:len_avr,:));
SE_PS_LSAA_2b_avr=mean(SE_PS_LSAA_2b(1:len_avr,:));
SE_PS_LSAA_4b_avr=mean(SE_PS_LSAA_4b(1:len_avr,:));

SW_rs_avr=mean(SW_rs(1:len_avr,:));
SW_random_avr=mean(SW_random(1:len_avr,:)); 

SW_Apga_TS=SW_apga.ts{1};
SW_Rdn_TS=SW_rdn.ts{1};
SW_Apga_TS_avr=mean(SW_Apga_TS(1:len_avr,:)); 
SW_Rdn_TS_avr=mean(SW_Rdn_TS(1:len_avr,:)); 

EE_DBF=EE_compute_DC(SE_DBF_avr,Nr,Nt);
EE_FTTD=EE_compute_FTTD_HBF(SE_FTTD_avr,Nr,Nt,N_rf,Nfdt,Nfdr);
EE_PS_LSAA_1b=EE_compute_PS_HBF(SE_PS_LSAA_1b_avr,Nr,Nt,N_rf,1);
EE_PS_LSAA_2b=EE_compute_PS_HBF(SE_PS_LSAA_2b_avr,Nr,Nt,N_rf,2);
EE_PS_LSAA_4b=EE_compute_PS_HBF(SE_PS_LSAA_4b_avr,Nr,Nt,N_rf,4);

EE_TTD_1b=EE_compute_TTD_HBF(SE_TTD_1b_avr,Nr,Nt,N_rf,Ndt,Ndr,1);
EE_TTD_2b=EE_compute_TTD_HBF(SE_TTD_2b_avr,Nr,Nt,N_rf,Ndt,Ndr,2);
EE_TTD_4b=EE_compute_TTD_HBF(SE_TTD_4b_avr,Nr,Nt,N_rf,Ndt,Ndr,4);

EE_SW_rs=EE_compute_SW_HBF(SW_rs_avr,Nr,Nt,N_rf);
SW_apga_EE=EE_compute_SW_HBF(SW_Apga_TS_avr,Nr,Nt,N_rf);
SW_Rdn_EE=EE_compute_SW_HBF(SW_Rdn_TS_avr,Nr,Nt,N_rf);



time_mark=datestr(now,'mmmm-dd');
% file_name=[saved_name,'_',time_mark];
file_name=['SE_BW1_LSAA',time_mark];
save(file_name,'-v7.3')
   
altstyles = {'m--o'; 'r-p';'g-h'; 'b-s'; 'c-v'; 'k-*'; 'r-.p';'g-.h'; 'b-.s'; 'c-.v'; 'k-.*'};
x=B_all/GHz;
% x=BSR_all;
% figure
% hold on
% h1=plot(x,SE_PS_pp_4b_avr,'b-o');
% h2=plot(x,SE_PS_pp_2b_avr,'b--o');
% h3=plot(x,SE_PS_pp_1b_avr,'b-.o');
% 
% h4=plot(x,SE_PS_LSAA_4b_avr,'r-s');
% h5=plot(x,SE_PS_LSAA_2b_avr,'r--s');
% h6=plot(x,SE_PS_LSAA_1b_avr,'r-.s');
% legend([h1,h2,h3,h4,h5,h6],'PS-HBF-pp:4bits','PS-HBF-pp:2bits','PS-HBF-pp:1bit','PS-HBF-LSAA:4bits','PS-HBF-LSAA:2bits','PS-HBF-LSAA:1bit');
% ylabel('Average SE (bits/s/Hz)')
% xlabel('BSR')
% % xticks(Nall)
% xlim([x(1) x(end)])
% grid on
% box on

figure
hold on
h0=plot(x,SE_DBF_avr,'k-*','LineWidth',1.5);
h1=plot(x,SE_FTTD_avr,'b-v','LineWidth',1.5);
h2=plot(x,SE_PS_LSAA_4b_avr,'m-o','LineWidth',1.5);
h3=plot(x,SE_PS_LSAA_2b_avr,'m--o','LineWidth',1.5);
h4=plot(x,SE_PS_LSAA_1b_avr,'m:o','LineWidth',1.5);

h5=plot(x,SE_TTD_4b_avr,'-s','LineWidth',1.5);
h6=plot(x,SE_TTD_2b_avr,'--s','LineWidth',1.5);
h7=plot(x,SE_TTD_1b_avr,':s','LineWidth',1.5);

h8=plot(x,SW_Apga_TS_avr,'-d','LineWidth',1.5);
h9=plot(x,SW_Rdn_TS_avr,'--d','LineWidth',1.5);
h10=plot(x,SW_rs_avr,':d','LineWidth',1.5);
h11=plot(x,SW_random_avr,'-x','LineWidth',1.5);


set(h5,'Color',[0,0.5,0])%牛油绿
set(h6,'Color',[0,0.5,0])%牛油绿
set(h7,'Color',[0,0.5,0])%牛油绿

set(h8,'Color',[0.6,0.2,0])%brown
set(h9,'Color',[0.6,0.2,0])%brown
set(h10,'Color',[0.6,0.2,0])%brown
set(h11,'Color',[0.4940 0.1840 0.5560])% purple

legend([h0,h1,h2,h3,h4,h5,h6,h7],'DBF','FTTD-HBF','PS-HBF:4bits','PS-HBF:2bits','PS-HBF:1bit','TTD-HBF:4bits','TTD-HBF:2bits','TTD-HBF:1bit');

ylabel('Average SE (bits/s/Hz)')
xlabel('Bandwidth (GHz)')
% xticks(Nall)
xlim([x(1) x(end)])
grid on
box on

axesNew = axes('position',get(gca,'position'),'visible','off');
legend(axesNew,[h8,h9,h10,h11],'SW-HBF:PGA-TS','SW-HBF:TS','SW-HBF:PGA-TbRs','SW-HBF:random','location','west')

figure
hold on
h0=plot(x,EE_DBF,'k-*','LineWidth',1.5);
h1=plot(x,EE_FTTD,'b-v','LineWidth',1.5);
h2=plot(x,EE_PS_LSAA_4b,'m-o','LineWidth',1.5);
h3=plot(x,EE_PS_LSAA_2b,'m--o','LineWidth',1.5);
h4=plot(x,EE_PS_LSAA_1b,'m-.o','LineWidth',1.5);

h5=plot(x,EE_TTD_4b,'-s','LineWidth',1.5);
h6=plot(x,EE_TTD_2b,'--s','LineWidth',1.5);
h7=plot(x,EE_TTD_1b,':s','LineWidth',1.5);

h8=plot(x,SW_apga_EE,'-d','LineWidth',1.5);
h9=plot(x,SW_Rdn_EE,'--d','LineWidth',1.5);
h10=plot(x,EE_SW_rs,':d','LineWidth',1.5);
% h11=plot(x,EE_SW_random);

set(h5,'Color',[0,0.5,0])%牛油绿
set(h6,'Color',[0,0.5,0])%牛油绿
set(h7,'Color',[0,0.5,0])%牛油绿

set(h8,'Color',[0.6,0.2,0])%brown
set(h9,'Color',[0.6,0.2,0])%brown
set(h10,'Color',[0.6,0.2,0])%brown
% set(h11,'Color',[0.4940 0.1840 0.5560])% purple

legend([h0,h1,h2,h3,h4,h5,h6,h7],'DBF','FTTD-HBF','PS-HBF:4bits','PS-HBF:2bits','PS-HBF:1bit','TTD-HBF:4bits','TTD-HBF:2bits','TTD-HBF:1bit');

ylabel('Average EE (Gbits/J)')
xlabel('Bandwidth')
% xticks(Nall)
xlim([x(1) x(end)])
grid on 
box on
axesNew = axes('position',get(gca,'position'),'visible','off');
legend(axesNew,[h8,h9,h10],'SW-HBF:PGA-TS','SW-HBF:TS','SW-HBF:PGA-TbRs','location','west')
figure
% xx=B_all/GHz;
plot(x,BSR_all,'r-o')
ylabel('BSR')
xlabel('Bandwidth (GHz)')
xticks(x)
xlim([x(1) x(end)])
grid on 
box on
ccc=1;

end
%% defined funcions
function EE_DC=EE_compute_DC(SE_DC,Nr,Nt)
factor2=1e-3;
P_LNA=60*factor2;
P_SP=19.5*factor2;
P_C=19.5*factor2;

P_M=22*factor2;
P_LO=4*factor2;
P_LPF=14*factor2;
P_BBamp=5*factor2;

P_ADC=88*factor2;
P_PS=40*factor2;
P_SW=5*factor2;
P_RF=P_M+P_LO+P_LPF+P_BBamp;

Consp_DC=(Nr+Nt)*(P_LNA+P_RF+2*P_ADC);

EE_DC=SE_DC/Consp_DC;
end

function EE_PS_HBF=EE_compute_PS_HBF(SE_PS_HBF,Nr,Nt,N_rf,bit_ps)
factor2=1e-3;
P_LNA=60*factor2;
P_SP=19.5*factor2;
P_C=19.5*factor2;

P_M=22*factor2;
P_LO=4*factor2;
P_LPF=14*factor2;
P_BBamp=5*factor2;

P_ADC=88*factor2;
switch bit_ps
    case 4
        P_PS=40*factor2;
    case 2
        P_PS=20*factor2; 
    case 1
        P_PS=10*factor2;     
end
P_RF=P_M+P_LO+P_LPF+P_BBamp;

part1=(Nt+Nr)*(P_LNA+N_rf*P_PS)+2*N_rf*(P_RF+2*P_ADC);
part2=(Nr+N_rf)*P_SP + (Nt+N_rf)*P_C;

Consp_PS=part1+part2;

EE_PS_HBF=SE_PS_HBF/Consp_PS;
end

function EE_SW_HBF=EE_compute_SW_HBF(SE_SW_HBF,Nr,Nt,N_rf)
factor2=1e-3;
P_LNA=60*factor2;
P_SP=19.5*factor2;
P_C=19.5*factor2;

P_M=22*factor2;
P_LO=4*factor2;
P_LPF=14*factor2;
P_BBamp=5*factor2;

P_ADC=88*factor2;
P_SW=5*factor2;
P_RF=P_M+P_LO+P_LPF+P_BBamp;

part1=(Nt+Nr)*(P_LNA+N_rf*P_SW)+2*N_rf*(P_RF+2*P_ADC);
part2=(Nr+N_rf)*P_SP + (Nt+N_rf)*P_C;

Consp_SW=part1+part2;

EE_SW_HBF=SE_SW_HBF./Consp_SW;
end

function EE_TTD_HBF=EE_compute_TTD_HBF(SE_TTD_HBF,Nr,Nt,N_rf,Ndt,Ndr,bit_ps)
factor2=1e-3;
P_LNA=60*factor2;
P_SP=19.5*factor2;
P_C=19.5*factor2;

P_M=22*factor2;
P_LO=4*factor2;
P_LPF=14*factor2;
P_BBamp=5*factor2;

P_ADC=88*factor2;
switch bit_ps
    case 4
        P_PS=40*factor2;
    case 2
        P_PS=20*factor2; 
    case 1
        P_PS=10*factor2;     
end
P_RF=P_M+P_LO+P_LPF+P_BBamp;

P_TTD=285*factor2;% ttd power consumption

part1=(Nt+Nr)*(P_LNA+N_rf*P_PS)+2*N_rf*(P_RF+2*P_ADC)+N_rf*(Ndt+Ndr)*P_TTD;
part2=(Nr+N_rf+N_rf*Ndt)*P_SP + (Nt+N_rf+N_rf*Ndr)*P_C;

Consp_TTD=part1+part2;

EE_TTD_HBF=SE_TTD_HBF/Consp_TTD;
end

function EE_TTD_HBF=EE_compute_FTTD_HBF(SE_FTTD_HBF,Nr,Nt,N_rf,Nfdt,Nfdr)
factor2=1e-3;
P_LNA=60*factor2;
P_SP=19.5*factor2;
P_C=19.5*factor2;

P_M=22*factor2;
P_LO=4*factor2;
P_LPF=14*factor2;
P_BBamp=5*factor2;

P_ADC=88*factor2;
P_SW=5*factor2;
P_RF=P_M+P_LO+P_LPF+P_BBamp;

P_FTTD=30*factor2;% ttd power consumption

part1=(Nt+Nr)*(P_LNA+P_SW)+2*N_rf*(P_RF+2*P_ADC)+N_rf*(Nfdt+Nfdr)*P_FTTD;
part2=N_rf*(Nfdt*P_SP + Nfdr*P_C);

Consp_FTTD=part1+part2;

EE_TTD_HBF=SE_FTTD_HBF/Consp_FTTD;
end