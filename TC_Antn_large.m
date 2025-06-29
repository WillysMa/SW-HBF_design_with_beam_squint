function TC=TC_Antn_large
clc;clear all;close all;
tic
rng(0)
%% variable: antennas
GHz=1e+9;%bandwidth=1GHz
Nc=64;% carrier number
D=Nc/4;% maximum delay
B=10*GHz;% bandwidth
fc=140*GHz;% carrier frequency
% Nt=256;
% Nr=Nt;
N_rf=4;
Ns=N_rf;
AtnDis=1/2;
Frac_Bw=B/fc;

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
PhaseSet_1bit=phase_Set{1};
PhaseSet_2bit=phase_Set{2};
PhaseSet_4bit=phase_Set{3};

Ndt=1;
Ndr=Ndt;
Nfdt=2;
Nfdr=Nfdt;

delta=0.1;
params_pga.acc=1e-4;
params_pga.iter_max=1000;%parameters for pga

Iter_max=200;
params_ts.iter_max=Iter_max;
params_ts.local_iter_max=10;%TS stop criteria
params_ts.Len_tl=Iter_max;% length of tabu list

params_TbRs.iter_max=Iter_max;
params_TbRs.local_iter_max=10;
params_TbRs.Len_tl=Iter_max;

N_all=2.^(3:8);
% N_all=10:70:500;
BSR_all=AtnDis*N_all*Frac_Bw/8;

len=length(N_all);
avr=100;


%% 
IniMatrix=zeros(avr,len);
SE_DBF=IniMatrix;

TC_ES_all=IniMatrix;
TC_FcPs_all=IniMatrix;
TC_DynPs_all=IniMatrix;
TC_FcTd_all=IniMatrix;
TC_DynFtd_all=IniMatrix;
TC_SW_all=IniMatrix;

obj_fun_pg=@compute_rate_v0;
obj_fun=@compute_rate_v0;        
parfor n=1:avr

    for id=1:len
        Nt=N_all(id);
        Nr=Nt;
    

        [Hcc,ch_params]=Channel_Gen_compact(Nt,Nr,L,AtnDis,B,ofdm_paras);
        
        [~,~,DBF,Proj_Comb]=DBF_water_filling(Hcc,Pb,Ns,noise_power);
 

        % [SE_ts,W_rf,tc_es]=SW_HBF_ES(Pb,Ns,N_rf,Hcc,noise_power,episilon,obj_fun);
        % TC_ES_all(n,id)=tc_es;

        %  [SE_PS_pp_1b(n,id),~]=HBF_pp_quantized(Hcc,Pb,Ns,N_rf,noise_power,episilon,phase_Set{1}); 
        %  [SE_PS_pp_2b(n,id),~]=HBF_pp_quantized(Hcc,Pb,Ns,N_rf,noise_power,episilon,phase_Set{2});   
        %  [SE_PS_pp_4b(n,id),~]=HBF_pp_quantized(Hcc,Pb,Ns,N_rf,noise_power,episilon,phase_Set{3});  
        % 
        %  obj_fun=@compute_rate_v0;        
        %  [SW_apga.ts{1}(n,id),~]=HBF_APGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,delta,obj_fun,obj_fun_pg);
        %  [SW_rdn.ts{1}(n,id),~]=HBF_TS(params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,obj_fun);      
        %  [SW_npga.ts{1}(n,id),~]=HBF_NPGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,obj_fun,obj_fun_pg);
        %  [SW_rs(n,id),~]=HBF_APGA_TbRs(params_pga,params_TbRs,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,delta,obj_fun,obj_fun_pg);  
        % 
        % obj_fun=@compute_rate_v1;
        %  [SW_apga.ts{2}(n,id),~]=HBF_APGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,delta,obj_fun,obj_fun_pg);
        %  [SW_rdn.ts{2}(n,id),~]=HBF_TS(params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,obj_fun);         
        %  [SW_npga.ts{2}(n,id),~]=HBF_NPGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,obj_fun,obj_fun_pg);
        % 
        %  [SW_random(n,id),~]=SW_HBF_random(Pb,Ns,N_rf,Hcc,Ini1,noise_power,episilon,obj_fun);
      
         % [SE_PS_LSAA_1b(n,id),~]=HBF_LSAA(Hcc,Pb,Ns,N_rf,noise_power,episilon,PhaseSet_1bit); 
%          [SE_PS_LSAA_2b(n,id),~]=HBF_LSAA(Hcc,Pb,Ns,N_rf,noise_power,episilon,PhaseSet_2bit);   
         [~,~,tc_fcps]=HBF_LSAA(Hcc,Pb,Ns,N_rf,noise_power,episilon,phase_Set{3});
         TC_FcPs_all(n,id)=tc_fcps;

         [~,~,tc_fctd]=HBF_DDP(Hcc,Pb,Ns,N_rf,Ndt,Ndr,noise_power,AtnDis,B,ch_params,ofdm_paras,PhaseSet_4bit);
         TC_FcTd_all(n,id)=tc_fctd;

         [~,~,tc_dynftd]=HBF_FTTD(Hcc,DBF,N_rf,Nfdt,Nfdr,noise_power,AtnDis,B,ofdm_paras);
         TC_DynFtd_all(n,id)=tc_dynftd;

         Flag_Identity_Ini=1;         
         [~,~,tc_swpp]=HBF_APGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,noise_power,episilon,delta,obj_fun,obj_fun_pg,Flag_Identity_Ini,Proj_Comb);
%          [SW_Rdn_TS(n,id),~]=HBF_TS(params_ts,Pb,Ns,N_rf,Hcc,noise_power,episilon,obj_fun,Flag_Identity_Ini,Proj_Comb);      
         % [SW_rs(n,id),~]=HBF_APGA_TbRs(params_pga,params_TbRs,Pb,Ns,N_rf,Hcc,noise_power,episilon,delta,obj_fun,obj_fun_pg,Flag_Identity_Ini,Proj_Comb); 
         TC_SW_all(n,id)=tc_swpp;

         iter_max=500;tol_obj=1e-4;
         [~,tc_dynps]=Dyn_PS_HBF(Hcc,N_rf,Ns,Pb,noise_power,PhaseSet_4bit,iter_max,tol_obj); 
         TC_DynPs_all(n,id)=tc_dynps;
        
    disp(['id=',num2str(id),', avr=',num2str(n)])
    % time_mark=datestr(now,'mmmm-dd');
    % file_name=['SE_Antn_',time_mark];
    % save(file_name,'-v7.3')
        ccc=1;
    end
end




TC_ES=mean(TC_ES_all);
TC_FcPs=mean(TC_FcPs_all);
TC_DynPs=mean(TC_DynPs_all);
TC_FcTd=mean(TC_FcTd_all);
TC_DynFtd=mean(TC_DynFtd_all);
TC_SW=mean(TC_SW_all);

TC=toc/3600;
time_mark=datestr(now,'mmmm-dd');
file_name=['TC_Antn_Large',time_mark];
save(file_name,'-v7.3')

altstyles = {'m--o'; 'r-p';'g-h'; 'b-s'; 'c-v'; 'k-*'; 'r-.p';'g-.h'; 'b-.s'; 'c-.v'; 'k-.*'};
x=N_all;
figure
hold on
h1=plot(x,TC_ES,altstyles{1},'LineWidth',1.5);
h2=plot(x,TC_FcPs,altstyles{2},'LineWidth',1.5);
h3=plot(x,TC_DynPs,altstyles{3},'LineWidth',1.5);
h4=plot(x,TC_FcTd,altstyles{4},'LineWidth',1.5);
h5=plot(x,TC_DynFtd,altstyles{5},'LineWidth',1.5);
h6=plot(x,TC_SW,altstyles{6},'LineWidth',1.5);

% set(h10,'Color',[0,0.5,0],'Marker','d','LineStyle','-')%牛油绿
% set(h11,'Color',[0.6,0.2,0],'Marker','^','LineStyle','-')%brown
% set(h0,'Color',[0.4940 0.1840 0.5560],'Marker','x','LineStyle','-')% purple

legend([h1,h2,h3,h4,h5,h6],'SW-HBF:ES','PS-HBF','Dyn-PS-HBF','FC-TTD-HBF','Dyn-FTTD-HBF','SW-HBF:PGA-TS');
ylabel('Time cost (s)')
xlabel('Number of antennas')
% xticks(Nall)
xlim([x(1) x(end)])
grid on
box on
ccc=1;
% SE_DBF_avr=mean(SE_DBF(1:len_avr,:));
% SE_PS_pp_1b_avr=mean(SE_PS_pp_1b(1:len_avr,:));
% SE_PS_pp_2b_avr=mean(SE_PS_pp_2b(1:len_avr,:));
% SE_PS_pp_4b_avr=mean(SE_PS_pp_4b(1:len_avr,:));
% 
% SW_rs_avr=mean(SW_rs(1:len_avr,:));
% SW_random_avr=mean(SW_random(1:len_avr,:)); 
% for i=1:2
%     SW_apga.ts_avr{i}=mean(SW_apga.ts{i}(1:len_avr,:)); 
%     SW_npga.ts_avr{i}=mean(SW_npga.ts{i}(1:len_avr,:)); 
%     SW_rdn.ts_avr{i}=mean(SW_rdn.ts{i}(1:len_avr,:)); 
% end

% SE_DBF_avr=SE_DBF(1:len_avr,1:x_num);
% SE_PS_pp_2b_avr=SE_PS_pp_2b(1:len_avr,1:x_num);
% SE_PS_pp_4b_avr=SE_PS_pp_4b(1:len_avr,1:x_num);
% 
% SW_rs_avr=SW_rs(1:len_avr,1:x_num);
% SW_random_avr=SW_random(1:len_avr,1:x_num); 
% for i=1:2
%     SW_apga.ts_avr{i}=SW_apga.ts{i}(1:len_avr,1:x_num); 
%     SW_npga.ts_avr{i}=SW_npga.ts{i}(1:len_avr,1:x_num); 
%     SW_rdn.ts_avr{i}=SW_rdn.ts{i}(1:len_avr,1:x_num); 
% end

% EE_DBF=EE_compute_DC(SE_DBF_avr,N_all,N_all);
% EE_PS_pp_1b=EE_compute_PS_HBF(SE_PS_pp_1b_avr,N_all,N_all,N_rf,1);
% EE_PS_pp_2b=EE_compute_PS_HBF(SE_PS_pp_2b_avr,N_all,N_all,N_rf,2);
% EE_PS_pp_4b=EE_compute_PS_HBF(SE_PS_pp_4b_avr,N_all,N_all,N_rf,4);
% EE_SW_rs=EE_compute_SW_HBF(SW_rs_avr,N_all,N_all,N_rf);
% 
% for i=1:2
%     SW_apga.EE{i}=EE_compute_SW_HBF(SW_apga.ts_avr{i},N_all,N_all,N_rf);
%     SW_npga.EE{i}=EE_compute_SW_HBF(SW_npga.ts_avr{i},N_all,N_all,N_rf);
%     SW_rdn.EE{i}=EE_compute_SW_HBF(SW_rdn.ts_avr{i},N_all,N_all,N_rf);
% end

% time_mark=datestr(now,'mmmm-dd');
% file_name=['SE_Antn_',time_mark];
%  save(file_name,'-v7.3')
   
% altstyles = {'m--o'; 'r-p';'g-h'; 'b-s'; 'c-v'; 'k-*'; 'r-.p';'g-.h'; 'b-.s'; 'c-.v'; 'k-.*'};
% % x=B_all/GHz;
% x=N_all;
% figure
% hold on
% % h0=plot(x,SE_DBF_avr);
% h1=plot(x,SE_PS_pp_4b_avr,'m-o');
% h2=plot(x,SE_PS_pp_2b_avr,'m--o');
% h3=plot(x,SE_PS_pp_1b_avr,'m-.o');
% 
% h4=plot(x,SW_apga.ts_avr{1},altstyles{2});
% h5=plot(x,SW_apga.ts_avr{2},altstyles{7});
% h6=plot(x,SW_npga.ts_avr{1},altstyles{3});
% h7=plot(x,SW_npga.ts_avr{2},altstyles{8});
% h8=plot(x,SW_rdn.ts_avr{1},altstyles{4});
% h9=plot(x,SW_rdn.ts_avr{2},altstyles{9});
% h10=plot(x,SW_rs_avr);
% h11=plot(x,SW_random_avr);
% 
% 
% set(h10,'Color',[0,0.5,0],'Marker','d','LineStyle','-')%牛油绿
% set(h11,'Color',[0.6,0.2,0],'Marker','^','LineStyle','-')%brown
% % set(h0,'Color',[0.4940 0.1840 0.5560],'Marker','x','LineStyle','-')% purple
% 
% legend([h1,h2,h3],'PS-HBF:4bits','PS-HBF:2bits','PS-HBF:1bit');
% ylabel('Average SE (bits/s/Hz)')
% xlabel('Number of antennas')
% % xticks(Nall)
% xlim([x(1) x(end)])
% grid on
% box on
% axesNew = axes('position',get(gca,'position'),'visible','off');
% legend(axesNew,[h4,h5,h6,h7,h8,h9,h10,h11],'SW-HBF:APGA-TS0','SW-HBF:APGA-TS1','SW-HBF:NPGA-TS0','SW-HBF:NPGA-TS1',...
%     'SW-HBF:TS0','SW-HBF:TS1 ','SW-HBF:APGA-TbRs','SW-HBF:random','location','west')
% 
% figure
% hold on
% % h0=plot(x,EE_DBF);
% h1=plot(x,EE_PS_pp_4b,'m-o');
% h2=plot(x,EE_PS_pp_2b,'m--o');
% h3=plot(x,EE_PS_pp_1b,'m-.o');
% 
% h4=plot(x,SW_apga.EE{1},altstyles{2});
% h5=plot(x,SW_apga.EE{2},altstyles{7});
% h6=plot(x,SW_npga.EE{1},altstyles{3});
% h7=plot(x,SW_npga.EE{2},altstyles{8});
% h8=plot(x,SW_rdn.EE{1},altstyles{4});
% h9=plot(x,SW_rdn.EE{2},altstyles{9});
% h10=plot(x,EE_SW_rs);
% 
% 
% set(h10,'Color',[0,0.5,0],'Marker','d','LineStyle','-')%牛油绿
% % set(h9,'Color',[0.6,0.2,0],'Marker','^','LineStyle','-')%brown
% % set(h0,'Color',[0.4940 0.1840 0.5560],'Marker','x','LineStyle','-')% purple
% 
% legend([h1,h2,h3],'PS-HBF:4bits','PS-HBF:2bits','PS-HBF:1bit');
% ylabel('Average EE (Gbits/J)')
% xlabel('Number of antennas')
% % xticks(Nall)
% xlim([x(1) x(end)])
% grid on 
% box on
% axesNew = axes('position',get(gca,'position'),'visible','off');
% legend(axesNew,[h4,h5,h6,h7,h8,h9,h10],'SW-HBF:APGA-TS0','SW-HBF:APGA-TS1','SW-HBF:NPGA-TS0','SW-HBF:NPGA-TS1',...
%     'SW-HBF:TS0','SW-HBF:TS1 ','SW-HBF:APGA-TbRs','location','west')
ccc=1;



