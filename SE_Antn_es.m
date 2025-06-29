clc;clear all;close all
%% variable: Nt
GHz=1e+9;%bandwidth=1GHz
Nc=64;% carrier number
D=Nc/4;% maximum delay
B=10*GHz;% bandwidth
fc=140*GHz;% carrier frequency
% Nt=4;
% Nr=Nt;
N_rf=2;
Ns=N_rf;
AtnDis=1/2;
Frac_Bw=B/fc;

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

params_ts.iter_max=200;
params_ts.local_iter_max=10;%TS stop criteria
params_ts.Len_tl=200;% length of tabu list

params_TbRs.iter_max=200;
params_TbRs.local_iter_max=10;
params_TbRs.Len_tl=params_TbRs.iter_max;

N_all=4:2:12;
BSR_coef=AtnDis*N_all*Frac_Bw/4;

len=length(N_all);
avr=50;

%% 
IniMatrix=zeros(avr,len);
SE_DBF=IniMatrix;

SE_PS_pp_qt=IniMatrix;
SE_PS_pp2_qt=IniMatrix;

for i=1:4
SW_pga.ts{i}=IniMatrix;
SW_pga.tc{i}=IniMatrix;
SW_rdn.ts{i}=IniMatrix;
SW_rdn.tc{i}=IniMatrix;
end

SW_rs=IniMatrix;
Rcd_rs=IniMatrix;

SW_es=IniMatrix;
Rcd_es=IniMatrix;

SW_pga_es=IniMatrix;
Rcd_pga_es=IniMatrix;

SW_random=IniMatrix;
Rcd_random=IniMatrix;
for n=1:avr

    for id=1:len 
        Nt=N_all(id);
        Nr=Nt;
        [Hcc,~]=Channel_Gen_compact(Nt,Nr,L,AtnDis,ofdm_paras);
        
        SE_DBF(n,id)=DBF_water_filling(Hcc,Pb,Ns,noise_power,episilon);
 
         [SE_PS_pp_qt(n,id),W_rf_ps]=HBF_pp_quantized(Hcc,Pb,Ns,N_rf,noise_power,episilon,phase_Set);   
%          [SE_PS_pp2_qt(n,id),W_rf_ps2]=HBF_pp2_quantized(Hcc,Pb,Ns,N_rf,noise_power,episilon,phase_Set);  

        tic
          t0=clock;
         obj_fun=@compute_rate_v0;
         obj_fun_pg=@compute_rate_v0;
         [SW_pga.ts{1}(n,id),W_rf0_pg]=SW_HBF_Tabu_PGA(params_pga,params_ts,Pb,Ns,N_rf,Hcc,noise_power,episilon,delta,obj_fun,obj_fun_pg);
         t0_1=clock;
         [SW_rdn.ts{1}(n,id),W_rf0]=SW_HBF_Tabu(params_ts,Pb,Ns,N_rf,Hcc,noise_power,episilon,obj_fun);


        t1=clock;
        obj_fun=@compute_rate_v1;
        obj_fun_pg=@compute_rate_v0;
         [SW_pga.ts{2}(n,id),W_rf1_pg]=SW_HBF_Tabu_PGA(params_pga,params_ts,Pb,Ns,N_rf,Hcc,noise_power,episilon,delta,obj_fun,obj_fun_pg);
         t1_1=clock;
         [SW_rdn.ts{2}(n,id),W_rf1]=SW_HBF_Tabu(params_ts,Pb,Ns,N_rf,Hcc,noise_power,episilon,obj_fun);

         t2=clock;

         obj_fun=@compute_rate_v2;
         obj_fun_pg=@compute_rate_v0;
         [SW_pga.ts{3}(n,id),W_rf2_pg]=SW_HBF_Tabu_PGA(params_pga,params_ts,Pb,Ns,N_rf,Hcc,noise_power,episilon,delta,obj_fun,obj_fun_pg);
          t2_1=clock;
         [SW_rdn.ts{3}(n,id),W_rf2]=SW_HBF_Tabu(params_ts,Pb,Ns,N_rf,Hcc,noise_power,episilon,obj_fun);

         t3=clock;
          obj_fun=@compute_rate_v3;
          obj_fun_pg=@compute_rate_v0;
         [SW_pga.ts{4}(n,id),W_rf3_pg]=SW_HBF_Tabu_PGA(params_pga,params_ts,Pb,Ns,N_rf,Hcc,noise_power,episilon,delta,obj_fun,obj_fun_pg);
         t3_1=clock;
         [SW_rdn.ts{4}(n,id),W_rf3]=SW_HBF_Tabu(params_ts,Pb,Ns,N_rf,Hcc,noise_power,episilon,obj_fun);

         t4=clock;

         obj_fun=@compute_rate_v0;
         obj_fun_pg=@compute_rate_v0;
         [SW_rs(n,id),W_rf_rs]=SW_HBF_TbRs(params_pga,params_TbRs,Pb,Ns,N_rf,Hcc,noise_power,episilon,delta,obj_fun,obj_fun_pg);

         t5=clock;         
         obj_fun=@compute_rate_v0;
         [SW_es(n,id),W_rf_es]=SW_HBF_ES(Pb,Ns,N_rf,Hcc,noise_power,episilon,obj_fun);
         t6=clock;    

         obj_fun=@compute_rate_v0;
         obj_fun_pg=@compute_rate_v0;
         [SW_pga_es(n,id),W_rf_pga_es]=SW_HBF_PGA_ES(params_pga,Pb,Ns,N_rf,Hcc,noise_power,episilon,delta,obj_fun,obj_fun_pg); 
         
         t7=clock;    
         [SW_random(n,id),W_rf_rdn]=SW_HBF_random(Pb,Ns,N_rf,Hcc,noise_power,episilon,obj_fun);
         t8=clock;
         
        SW_pga.tc{1}(n,id)=etime(t0_1,t0);
        SW_rdn.tc{1}(n,id)=etime(t1,t0_1);

        SW_pga.tc{2}(n,id)=etime(t1_1,t1);
        SW_rdn.tc{2}(n,id)=etime(t2,t1_1);

        SW_pga.tc{3}(n,id)=etime(t2_1,t2);
        SW_rdn.tc{3}(n,id)=etime(t3,t2_1);

        SW_pga.tc{4}(n,id)=etime(t3_1,t3);
        SW_rdn.tc{4}(n,id)=etime(t4,t3_1);

        Rcd_rs(n,id)=etime(t5,t4);      
        Rcd_es(n,id)=etime(t6,t5);       
        Rcd_pga_es(n,id)=etime(t7,t6);
        Rcd_random(n,id)=etime(t8,t7);

    disp(['id=',num2str(id),', avr=',num2str(n)])
    time_mark=datestr(now,'mmmm-dd');
    file_name=['SE_Antn_es_',time_mark];
    save(file_name,'-v7.3')
        ccc=1;
    end
end

len=20;
SE_DBF_avr=mean(SE_DBF(1:len,:));
SE_PS_pp_qt_avr=mean(SE_PS_pp_qt(1:len,:));
% SE_PS_pp2_qt_avr=mean(SE_PS_pp2_qt(1:len,:));
for i=1:4
    SW_pga.ts_avr{i}=mean(SW_pga.ts{i}(1:len,:));
    SW_pga.tc_avr{i}=mean(SW_pga.tc{i}(1:len,:));
    SW_rdn.ts_avr{i}=mean(SW_rdn.ts{i}(1:len,:));
    SW_rdn.tc_avr{i}=mean(SW_rdn.tc{i}(1:len,:));
end
SW_rs_avr=mean(SW_rs(1:len,:));
Rcd_rs_avr=mean(Rcd_rs(1:len,:));

SW_es_avr=mean(SW_es(1:len,:));
Rcd_es_avr=mean(Rcd_es(1:len,:));

SW_pga_es_avr=mean(SW_pga_es(1:len,:));
Rcd_pga_es_avr=mean(Rcd_pga_es(1:len,:));

SW_random_avr=mean(SW_random(1:len,:));
Rcd_random_avr=mean(Rcd_random(1:len,:));
% time_mark=datestr(now,'mmmm-dd');
% file_name=['SE_Antn_es_',time_mark];
%  save(file_name,'-v7.3')
    
altstyles = {'r-o'; 'b-p';'b-.p'; 'g-s'; 'k-v'; 'm-*';  'g--s'; 'k--v'; 'm--*'; 'g-';'b-'; 'mv-'};
% x=B_all/GHz;
x=N_all;
figure
hold on
h1=plot(x,SE_DBF_avr,altstyles{1});
h2=plot(x,SE_PS_pp_qt_avr,altstyles{2});
h3=plot(x,SW_pga_es_avr,altstyles{3});

h4=plot(x,SW_pga.ts_avr{1},altstyles{4});
h5=plot(x,SW_pga.ts_avr{2},altstyles{5});
h6=plot(x,SW_pga.ts_avr{3},altstyles{6});
h7=plot(x,SW_pga.ts_avr{4});
h8=plot(x,SW_rs_avr);

h9=plot(x,SW_rdn.ts_avr{1},altstyles{7});
h10=plot(x,SW_rdn.ts_avr{2},altstyles{8});
h11=plot(x,SW_rdn.ts_avr{3},altstyles{9});
h12=plot(x,SW_rdn.ts_avr{4});

h13=plot(x,SW_es_avr,'r-.o');
h14=plot(x,SW_random_avr,'m-.>');

set(h7,'Color',[0,0.5,0],'Marker','d','LineStyle','-')%牛油绿
set(h8,'Color',[0.6,0.2,0],'Marker','^','LineStyle','-')%咖啡色
set(h12,'Color',[0,0.5,0],'Marker','d','LineStyle','--')%牛油绿

legend([h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,h13,h14],'Optimal DBF','PS-HBF-pp1','SW-HBF:PGA-ES','SW-HBF:PGA-TS0','SW-HBF:PGA-TS1','SW-HBF:PGA-TS2','SW-HBF:PGA-TS3','SW-HBF:PGA-TbRs',...
 'SW-HBF:TS0','SW-HBF:TS1','SW-HBF:TS2','SW-HBF:TS3','SW-HBF:ES','SW-HBF:random' )
ylabel('Average SE (bits/s/Hz)')
xlabel('Number of antannas $N_t=N_r$','interpreter','latex')
grid on 

figure
hold on
h1=semilogy(x,SW_pga.tc_avr{1},altstyles{4});
h2=semilogy(x,SW_pga.tc_avr{2},altstyles{5});
h3=semilogy(x,SW_pga.tc_avr{3},altstyles{6});
h4=semilogy(x,SW_pga.tc_avr{4});
h5=semilogy(x,Rcd_rs_avr);

h6=semilogy(x,SW_rdn.tc_avr{1},altstyles{7});
h7=semilogy(x,SW_rdn.tc_avr{2},altstyles{8});
h8=semilogy(x,SW_rdn.tc_avr{3},altstyles{9});
h9=semilogy(x,SW_rdn.tc_avr{4});

h10=semilogy(x,Rcd_es_avr,'r-.o');

h11=semilogy(x,Rcd_pga_es_avr,'r-.s');
h12=semilogy(x,Rcd_random_avr,'k-.>');
set(h4,'Color',[0,0.5,0],'Marker','d','LineStyle','-','LineWidth',1)%牛油绿
set(h5,'Color',[0.6,0.2,0],'Marker','^','LineStyle','-','LineWidth',1)%咖啡色
set(h9,'Color',[0,0.5,0],'Marker','d','LineStyle','--','LineWidth',1)%牛油绿

legend([h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12],'SW-HBF:PGA-TS0','SW-HBF:PGA-TS1','SW-HBF:PGA-TS2','SW-HBF:PGA-TS3','SW-HBF:PGA-TbRs',...
   'SW-HBF:TS0','SW-HBF:TS1','SW-HBF:TS2','SW-HBF:TS3', 'SW-HBF:ES','SW-HBF:PGA-ES','SW-HBF:random')

ylabel('Running time')
xlabel('Number of antannas $N_t=N_r$','interpreter','latex')
grid on 

% figure
% plot(B_all/GHz,BSR_all,'r-o')
% ylabel('BSR')
% xlabel('Bandwidth(GHz)')
% grid on 
ccc=1;


