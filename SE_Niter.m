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
AtnDis=1/2;

BSR_coef=AtnDis*Nt/fc/4;
ofdm_paras.B=B;
ofdm_paras.Nc=Nc;% carrier number
ofdm_paras.D=D;% maximum delay
ofdm_paras.fc=fc;% carrier frequency
L=5;

episilon=1e-4;

Pb_dBm=30;
Pb=10^(Pb_dBm/10);%mW

SNRx_dB=10;
SNRx=10.^(SNRx_dB/10);
noise_power=Pb./SNRx;% noise power


Phase_Bit=4;
phase_bit_num=length(Phase_Bit);
Nq_set=2.^Phase_Bit;
phase_Set=(-pi:2*pi/Nq_set:pi);

delta=0.1;
params_pga.acc=1e-4;
params_pga.iter_max=300;%parameters for pga




Niter_all=5:5:20;
lenN=length(Niter_all);

avr=20;
%%
IniMatrix=zeros(avr,lenN);
SE_DBF=zeros(1,avr);

SE_PS_pp_qt=zeros(1,avr);
SE_PS_pp2_qt=IniMatrix;


for i=1:3
    SW_rs.iter_max{i}=IniMatrix;
    SW_apga_ts.iter_max{i}=IniMatrix;
    SW_npga_ts.iter_max{i}=IniMatrix;
    SW_rdn_ts.iter_max{i}=IniMatrix;
end


obj_fun_pg=@compute_rate_v0;
obj_fun=@compute_rate_v0;
for n=1:avr
[Hcc,~]=Channel_Gen_compact(Nt,Nr,L,AtnDis,ofdm_paras);
%          [SE_DBF(n),DBF]=DBF_water_filling(Hcc,Pb,Ns,noise_power);
%          [SE_PS_pp_qt(n),W_rf_ps]=HBF_pp_quantized(Hcc,Pb,Ns,N_rf,noise_power,episilon,phase_Set); 
%          [SE_PS_pp2_qt(n,id),W_rf_ps2]=HBF_pp2_quantized(Hcc,Pb,Ns,N_rf,noise_power,episilon,phase_Set);  
    for id=1:lenN
        Niter=Niter_all(id);
        params_ts.local_iter_max=Niter;%TS stop criteria 
        params_TbRs.local_iter_max=Niter;
        
        Iter_max=100;
        params_ts.iter_max=Iter_max;      
        params_ts.Len_tl=Iter_max;% length of tabu list

        params_TbRs.iter_max=Iter_max;       
        params_TbRs.Len_tl=Iter_max;
    
         [SW_rdn_ts.iter_max{1}(n,id),~]=HBF_TS(params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,obj_fun);           
         [SW_apga_ts.iter_max{1}(n,id),~]=HBF_APGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,delta,obj_fun,obj_fun_pg);     
         [SW_npga_ts.iter_max{1}(n,id),~]=HBF_NPGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,obj_fun,obj_fun_pg);
         [SW_rs.iter_max{1}(n,id),~]=HBF_APGA_TbRs(params_pga,params_TbRs,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,delta,obj_fun,obj_fun_pg);
    

        Iter_max=200;
        params_ts.iter_max=Iter_max;      
        params_ts.Len_tl=Iter_max;% length of tabu list

        params_TbRs.iter_max=Iter_max;       
        params_TbRs.Len_tl=Iter_max;
        
         [SW_rdn_ts.iter_max{2}(n,id),~]=HBF_TS(params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,obj_fun);           
         [SW_apga_ts.iter_max{2}(n,id),~]=HBF_APGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,delta,obj_fun,obj_fun_pg);     
         [SW_npga_ts.iter_max{2}(n,id),~]=HBF_NPGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,obj_fun,obj_fun_pg);
         [SW_rs.iter_max{2}(n,id),~]=HBF_APGA_TbRs(params_pga,params_TbRs,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,delta,obj_fun,obj_fun_pg);
         
         
        Iter_max=300;
        params_ts.iter_max=Iter_max;      
        params_ts.Len_tl=Iter_max;% length of tabu list

        params_TbRs.iter_max=Iter_max;       
        params_TbRs.Len_tl=Iter_max;
        
         [SW_rdn_ts.iter_max{3}(n,id),~]=HBF_TS(params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,obj_fun);           
         [SW_apga_ts.iter_max{3}(n,id),~]=HBF_APGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,delta,obj_fun,obj_fun_pg);     
         [SW_npga_ts.iter_max{3}(n,id),~]=HBF_NPGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,obj_fun,obj_fun_pg);
         [SW_rs.iter_max{3}(n,id),~]=HBF_APGA_TbRs(params_pga,params_TbRs,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,delta,obj_fun,obj_fun_pg);
         
         
    disp(['id=',num2str(id),', avr=',num2str(n)])
    time_mark=datestr(now,'mmmm-dd');
    file_name=['SE_Niter_',time_mark];
    save(file_name,'-v7.3')
        ccc=1;
    end
    
end

len_avr=8;

for i=1:3
    SW_rs.avr{i}=mean(SW_rs.iter_max{i}(1:len_avr,:));
    SW_apga_ts.avr{i}=mean(SW_apga_ts.iter_max{i}(1:len_avr,:));
    SW_npga_ts.avr{i}=mean(SW_npga_ts.iter_max{i}(1:len_avr,:));
    SW_rdn_ts.avr{i}=mean(SW_rdn_ts.iter_max{i}(1:len_avr,:));
end

time_mark=datestr(now,'mmmm-dd');
file_name=['SE_Niter_',time_mark];
 save(file_name,'-v7.3')

altstyles = {'r-p';'g-h'; 'b-s'; 'k-o'; 'r-.p';'g-.h'; 'b-.s'; 'k-.o'; 'r--p';'g--h'; 'b--s'; 'k--o' };

x=Niter_all;
figure
hold on
h1=plot(x,SW_apga_ts.avr{1},altstyles{1});

h2=plot(x,SW_npga_ts.avr{1},altstyles{2});
h3=plot(x,SW_rdn_ts.avr{1},altstyles{3});
h4=plot(x,SW_rs.avr{1},altstyles{4});

h5=plot(x,SW_apga_ts.avr{2},altstyles{5});
h6=plot(x,SW_npga_ts.avr{2},altstyles{6});
h7=plot(x,SW_rdn_ts.avr{2},altstyles{7});
h8=plot(x,SW_rs.avr{2},altstyles{8});

h9=plot(x,SW_apga_ts.avr{3},altstyles{9});
h10=plot(x,SW_npga_ts.avr{3},altstyles{10});
h11=plot(x,SW_rdn_ts.avr{3},altstyles{11});
h12=plot(x,SW_rs.avr{3},altstyles{12});


% set(h7,'Color',[0,0.5,0],'Marker','d','LineStyle','-')%牛油绿
% set(h8,'Color',[0.6,0.2,0],'Marker','^','LineStyle','-')%brown
% set(h9,'Color',[0.4940 0.1840 0.5560],'Marker','x','LineStyle','-')% purple


legend([h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12],'APGA-TS0 $I=100$','NPGA-TS0 $I=100$','TS0 $I=100$','APGA-TbRs $I=100$',...
'APGA-TS0 $I=200$','NPGA-TS0 $I=200$','TS0 $I=200$','APGA-TbRs $I=200$','APGA-TS0 $3=100$','NPGA-TS0 $I=300$','TS0 $I=300$','APGA-TbRs $I=300$','interpreter','latex')
ylabel('Average SE (bits/s/Hz)')
xlabel('Niter')
grid on 
box on

cccc=1;