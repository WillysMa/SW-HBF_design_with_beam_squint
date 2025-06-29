clc;clear all;close all
%% variable: Bandwidth
GHz=1e+9;%bandwidth=1GHz
Nc=64;% carrier number
D=Nc/4;% maximum delay
B=10*GHz;% bandwidth
fc=140*GHz;% carrier frequency
Nt=256;
Nr=Nt;
N_rf=4;
Ns=N_rf;
% AtnDis=1/4;% ------------------------spatial oversampling------------------------------------
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


Phase_Bit=[1 2 4];
phase_bit_num=length(Phase_Bit);
Nq_set=2.^Phase_Bit;
for lb=1:phase_bit_num   
phase_Set{lb}=(-pi:2*pi/Nq_set(lb):pi);
end

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


AtnDis_all=2.^(-1*(5:-1:0));
BSR_all=AtnDis_all*Nt/8*B/fc;
len=length(AtnDis_all);

avr=100;

%% 
IniMatrix=zeros(avr,len);
SE_DBF=IniMatrix;

SE_PS_pp_4b=IniMatrix;
SE_PS_pp_2b=IniMatrix;
SE_PS_pp_1b=IniMatrix;

SW_rs=IniMatrix;
SW_random=IniMatrix;  
for i=1:2
    SW_apga.ts{i}=IniMatrix;
    SW_npga.ts{i}=IniMatrix;
    SW_rdn.ts{i}=IniMatrix;
end

for c=1:Nc
    Ini1.proj_comb{c}=eye(Nr);
end 
obj_fun_pg=@compute_rate_v0;        

for n=1:avr

    for id=1:len   
        AtnDis=AtnDis_all(id);
        [Hcc,~]=Channel_Gen_compact(Nt,Nr,L,AtnDis,ofdm_paras);
        
        [SE_DBF(n,id),DBF]=DBF_water_filling(Hcc,Pb,Ns,noise_power);
         [SE_PS_pp_1b(n,id),~]=HBF_pp_quantized(Hcc,Pb,Ns,N_rf,noise_power,episilon,phase_Set{1}); 
         [SE_PS_pp_2b(n,id),~]=HBF_pp_quantized(Hcc,Pb,Ns,N_rf,noise_power,episilon,phase_Set{2});   
         [SE_PS_pp_4b(n,id),~]=HBF_pp_quantized(Hcc,Pb,Ns,N_rf,noise_power,episilon,phase_Set{3});  

         obj_fun=@compute_rate_v0;        
         [SW_apga.ts{1}(n,id),~]=HBF_APGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,delta,obj_fun,obj_fun_pg);
         [SW_rdn.ts{1}(n,id),~]=HBF_TS(params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,obj_fun);      
%          [SW_npga.ts{1}(n,id),~]=HBF_NPGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,obj_fun,obj_fun_pg);

%         obj_fun=@compute_rate_v1;
%          [SW_apga.ts{2}(n,id),~]=HBF_APGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,delta,obj_fun,obj_fun_pg);
%          [SW_rdn.ts{2}(n,id),~]=HBF_TS(params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,obj_fun);         
%          [SW_npga.ts{2}(n,id),~]=HBF_NPGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,obj_fun,obj_fun_pg);
        
         obj_fun=@compute_rate_v0;
         [SW_rs(n,id),~]=HBF_APGA_TbRs(params_pga,params_TbRs,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,delta,obj_fun,obj_fun_pg);  
         [SW_random(n,id),~]=SW_HBF_random(Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,obj_fun);



    disp(['id=',num2str(id),', avr=',num2str(n)])
    time_mark=datestr(now,'mmmm-dd');
    file_name=['SE_AtDs_',time_mark];
    save(file_name,'-v7.3')
        ccc=1;
    end
end

len_avr=25;

SE_DBF_avr=mean(SE_DBF(1:len_avr,:));
SE_PS_pp_1b_avr=mean(SE_PS_pp_1b(1:len_avr,:));
SE_PS_pp_2b_avr=mean(SE_PS_pp_2b(1:len_avr,:));
SE_PS_pp_4b_avr=mean(SE_PS_pp_4b(1:len_avr,:));

SE_PS_LSAA_1b_avr=mean(SE_PS_pp_1b(1:len_avr,:));
SE_PS_LSAA_2b_avr=mean(SE_PS_pp_2b(1:len_avr,:));
SE_PS_LSAA_4b_avr=mean(SE_PS_pp_4b(1:len_avr,:));

SW_rs_avr=mean(SW_rs(1:len_avr,:));
SW_random_avr=mean(SW_random(1:len_avr,:)); 
for i=1:2
    SW_apga.ts_avr{i}=mean(SW_apga.ts{i}(1:len_avr,:)); 
    SW_npga.ts_avr{i}=mean(SW_npga.ts{i}(1:len_avr,:)); 
    SW_rdn.ts_avr{i}=mean(SW_rdn.ts{i}(1:len_avr,:)); 
end

EE_DBF=EE_compute_DC(SE_DBF_avr,Nr,Nt);
EE_PS_pp_1b=EE_compute_PS_HBF(SE_PS_pp_1b_avr,Nr,Nt,N_rf,1);
EE_PS_pp_2b=EE_compute_PS_HBF(SE_PS_pp_2b_avr,Nr,Nt,N_rf,2);
EE_PS_pp_4b=EE_compute_PS_HBF(SE_PS_pp_4b_avr,Nr,Nt,N_rf,4);
EE_SW_rs=EE_compute_SW_HBF(SW_rs_avr,Nr,Nt,N_rf);

for i=1:2
    SW_apga.EE{i}=EE_compute_SW_HBF(SW_apga.ts_avr{i},Nr,Nt,N_rf);
    SW_npga.EE{i}=EE_compute_SW_HBF(SW_npga.ts_avr{i},Nr,Nt,N_rf);
    SW_rdn.EE{i}=EE_compute_SW_HBF(SW_rdn.ts_avr{i},Nr,Nt,N_rf);
end

time_mark=datestr(now,'mmmm-dd');
file_name=['SE_AtDs_',time_mark];
 save(file_name,'-v7.3')
   
altstyles = {'m--o'; 'r-p';'g-h'; 'b-s'; 'c-v'; 'k-*'; 'r-.p';'g-.h'; 'b-.s'; 'c-.v'; 'k-.*'};
% x=B_all/GHz;
x=AtnDis_all;

figure
hold on
h0=plot(x,SE_DBF_avr);
h1=plot(x,SE_PS_pp_4b_avr,'m-o');
h2=plot(x,SE_PS_pp_2b_avr,'m--o');
h3=plot(x,SE_PS_pp_1b_avr,'m-.o');

h4=plot(x,SW_apga.ts_avr{1},altstyles{2});
% h5=plot(x,SW_apga.ts_avr{2},altstyles{7});
% h6=plot(x,SW_npga.ts_avr{1},altstyles{3});
% h7=plot(x,SW_npga.ts_avr{2},altstyles{8});
h8=plot(x,SW_rdn.ts_avr{1},altstyles{4});
% h9=plot(x,SW_rdn.ts_avr{2},altstyles{9});
h10=plot(x,SW_rs_avr);
h11=plot(x,SW_random_avr);



set(h10,'Color',[0,0.5,0],'Marker','d','LineStyle','-')%牛油绿
set(h11,'Color',[0.6,0.2,0],'Marker','^','LineStyle','-')%brown
set(h0,'Color',[0.4940 0.1840 0.5560],'Marker','x','LineStyle','-')% purple

legend([h0,h1,h2,h3],'DBF','PS-HBF:4bits','PS-HBF:2bits','PS-HBF:1bit');

ylabel('Average SE (bits/s/Hz)')
xlabel('Normalized antenna spacing')
xticks(x)
xlim([x(1) x(end)])
grid on
box on

axesNew = axes('position',get(gca,'position'),'visible','off');
legend(axesNew,[h4,h8,h10,h11],'SW-HBF:PGA-TS','SW-HBF:TS','SW-HBF:PGA-TbRs','SW-HBF:random','location','west')

 
figure
hold on
h0=plot(x,EE_DBF);
h1=plot(x,EE_PS_pp_4b,'m-o');
h2=plot(x,EE_PS_pp_2b,'m--o');
h3=plot(x,EE_PS_pp_1b,'m-.o');

h4=plot(x,SW_apga.EE{1},altstyles{2});
% h5=plot(x,SW_apga.EE{2},altstyles{7});
% h6=plot(x,SW_npga.EE{1},altstyles{3});
% h7=plot(x,SW_npga.EE{2},altstyles{8});
h8=plot(x,SW_rdn.EE{1},altstyles{4});
% h9=plot(x,SW_rdn.EE{2},altstyles{9});
h10=plot(x,EE_SW_rs);


set(h10,'Color',[0,0.5,0],'Marker','d','LineStyle','-')%牛油绿
% set(h11,'Color',[0.6,0.2,0],'Marker','^','LineStyle','-')%brown
set(h0,'Color',[0.4940 0.1840 0.5560],'Marker','x','LineStyle','-')% purple


legend([h0,h1,h2,h3],'DBF','PS-HBF:4bits','PS-HBF:2bits','PS-HBF:1bit');
ylabel('Average EE (Gbits/J)')
xlabel('Normalized antenna spacing')
xticks(x)
xlim([x(1) x(end)])
grid on 
box on

axesNew = axes('position',get(gca,'position'),'visible','off');
legend(axesNew,[h4,h8,h10],'SW-HBF:PGA-TS','SW-HBF:TS','SW-HBF:PGA-TbRs','location','west')

figure
plot(x,BSR_all,'r-o')
ylabel('BSR')
xlabel('Bandwidth')
xticks(x)
xlim([x(1) x(end)])
grid on 
box on
ccc=1;


%% defined funcions
function EE_DC=EE_compute_DC(SE_DC,Nr,Nt)
factor2=1e-3;
P_LNA=39*factor2;
P_SP=19.5*factor2;
P_C=19.5*factor2;

P_M=19*factor2;
P_LO=5*factor2;
P_LPF=14*factor2;
P_BBamp=5*factor2;

P_ADC=240*factor2;
P_PS=30*factor2;
P_SW=5*factor2;
P_RF=P_M+P_LO+P_LPF+P_BBamp;

Consp_DC=(Nr+Nt)*(P_LNA+P_RF+2*P_ADC);

EE_DC=SE_DC/Consp_DC;
end

function EE_PS_HBF=EE_compute_PS_HBF(SE_PS_HBF,Nr,Nt,N_rf,bit_ps)
factor2=1e-3;
P_LNA=39*factor2;
P_SP=19.5*factor2;
P_C=19.5*factor2;

P_M=19*factor2;
P_LO=5*factor2;
P_LPF=14*factor2;
P_BBamp=5*factor2;

P_ADC=240*factor2;
switch bit_ps
    case 4
        P_PS=40*factor2;
    case 2
        P_PS=20*factor2; 
    case 1
        P_PS=10*factor2;     
end
P_RF=P_M+P_LO+P_LPF+P_BBamp;


Consp_PS=(Nt+Nr)*(P_LNA+P_SP+N_rf*P_PS)+2*N_rf*(P_RF+P_C+2*P_ADC);

EE_PS_HBF=SE_PS_HBF/Consp_PS;
end

function EE_SW_HBF=EE_compute_SW_HBF(SE_SW_HBF,Nr,Nt,N_rf)
factor2=1e-3;
P_LNA=39*factor2;
P_SP=19.5*factor2;
P_C=19.5*factor2;

P_M=19*factor2;
P_LO=5*factor2;
P_LPF=14*factor2;
P_BBamp=5*factor2;

P_ADC=240*factor2;
P_PS=30*factor2;
P_SW=5*factor2;
P_RF=P_M+P_LO+P_LPF+P_BBamp;

Consp_SW=(Nt+Nr)*(P_LNA+P_SP+N_rf*P_SW)+2*N_rf*(P_RF+P_C+2*P_ADC);

EE_SW_HBF=SE_SW_HBF./Consp_SW;
end