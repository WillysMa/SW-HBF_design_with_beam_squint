clc;clear all;close all
%% variable: Bandwidth
GHz=1e+9;%bandwidth=1GHz
Nc=64;% carrier number
D=Nc/4;% maximum delay
% B=10*GHz;% bandwidth
fc=140*GHz;% carrier frequency
Nt=256;
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


BSR_all=0.1:0.3:2.2;
B_all=BSR_all/BSR_coef;
len=length(B_all);

avr=100;
SE_DBF_all=[];
SE_PS_pp_1b_all = [];
SE_PS_pp_2b_all = [];
SE_PS_pp_4b_all = [];
SW_rs_all = [];
SW_random_all = [];
SW_PGA_TS_all = [];
SW_PGA_TbRS_all = [];


EE_DBF_all=[];
EE_PS_pp_1b_all=[];
EE_PS_pp_2b_all=[];
EE_PS_pp_4b_all=[];
EE_SW_rs_all=[];
EE_PGA_TS_all = [];
EE_PGA_TbRS_all = [];

file_list={'SE_BW_July-02.mat', 'SE_BW_July-29.mat','SE_BW_August-18.mat'};%'SE_BW_June-13.mat'
for n=1:length(file_list)
    disp(file_list{n})
    load(file_list{n})
    SE_DBF_all=[SE_DBF_all;SE_DBF_avr];
    SE_PS_pp_1b_all=[SE_PS_pp_1b_all;SE_PS_pp_1b_avr];
    SE_PS_pp_2b_all=[SE_PS_pp_2b_all;SE_PS_pp_2b_avr];
    SE_PS_pp_4b_all=[SE_PS_pp_4b_all;SE_PS_pp_4b_avr];

    SW_rs_all=[SW_rs_all;SW_rs_avr];
    SW_random_all=[SW_random_all;SW_random_avr];

    SW_PGA_TS_all = [SW_PGA_TS_all; SW_apga.ts_avr{1}];
    SW_PGA_TbRS_all = [SW_PGA_TbRS_all;SW_rdn.ts_avr{1}];


    EE_DBF_all=[EE_DBF_all;EE_DBF];
    EE_PS_pp_1b_all=[EE_PS_pp_1b_all;EE_PS_pp_1b];
    EE_PS_pp_2b_all=[EE_PS_pp_2b_all;EE_PS_pp_2b];
    EE_PS_pp_4b_all=[EE_PS_pp_4b_all;EE_PS_pp_4b];
    EE_SW_rs_all=[EE_SW_rs_all;EE_SW_rs];

    EE_PGA_TS_all = [EE_PGA_TS_all;SW_apga.EE{1}];
    EE_PGA_TbRS_all = [EE_PGA_TbRS_all;SW_rdn.EE{1}];

end

SE_DBF_avr=mean(SE_DBF_all);
SE_PS_pp_1b_avr=mean(SE_PS_pp_1b_all);
SE_PS_pp_2b_avr=mean(SE_PS_pp_2b_all);
SE_PS_pp_4b_avr=mean(SE_PS_pp_4b_all);

SW_PGA_TS_avr=mean(SW_PGA_TS_all);
SW_PGA_TbRS_avr=mean(SW_PGA_TbRS_all);

SW_rs_avr=mean(SW_rs_all);
SW_random_avr=mean(SW_random_all); 

EE_DBF_avr=mean(EE_DBF_all);
EE_PS_pp_1b_avr=mean(EE_PS_pp_1b_all);
EE_PS_pp_2b_avr=mean(EE_PS_pp_2b_all);
EE_PS_pp_4b_avr=mean(EE_PS_pp_4b_all);

EE_PGA_TS_avr=mean(EE_PGA_TS_all);
EE_PGA_TbRS_avr=mean(EE_PGA_TbRS_all);
EE_SW_rs_avr=mean(EE_SW_rs_all);

altstyles = {'m--o'; 'r-p';'g-h'; 'b-s'; 'c-v'; 'k-*'; 'r-.p';'g-.h'; 'b-.s'; 'c-.v'; 'k-.*'};
x=BSR_all;

figure
hold on
h0=plot(x,SE_DBF_avr);
h1=plot(x,SE_PS_pp_4b_avr,'m-o');
h2=plot(x,SE_PS_pp_2b_avr,'m--o');
h3=plot(x,SE_PS_pp_1b_avr,'m-.o');

h4=plot(x,SW_PGA_TS_avr,altstyles{2});
% h5=plot(x,SW_apga.ts_avr{2},altstyles{7});
% h6=plot(x,SW_npga.ts_avr{1},altstyles{3});
% h7=plot(x,SW_npga.ts_avr{2},altstyles{8});
h8=plot(x,SW_PGA_TbRS_avr,altstyles{4});
% h9=plot(x,SW_rdn.ts_avr{2},altstyles{9});
h10=plot(x,SW_rs_avr);
h11=plot(x,SW_random_avr);



set(h10,'Color',[0,0.5,0],'Marker','d','LineStyle','-')%牛油绿
set(h11,'Color',[0.6,0.2,0],'Marker','^','LineStyle','-')%brown
set(h0,'Color',[0.4940 0.1840 0.5560],'Marker','x','LineStyle','-')% purple

legend([h0,h1,h2,h3],'DBF','PS-HBF:4bits','PS-HBF:2bits','PS-HBF:1bit');

ylabel('Average SE (bits/s/Hz)')
xlabel('BSR')
xticks(x)
xlim([x(1) x(end)])
grid on
box on

axesNew = axes('position',get(gca,'position'),'visible','off');
legend(axesNew,[h4,h8,h10,h11],'SW-HBF:PGA-TS','SW-HBF:TS','SW-HBF:PGA-TbRs','SW-HBF:random','location','west')

 
figure
hold on
h0=plot(x,EE_DBF_avr);
h1=plot(x,EE_PS_pp_4b_avr,'m-o');
h2=plot(x,EE_PS_pp_2b_avr,'m--o');
h3=plot(x,EE_PS_pp_1b_avr,'m-.o');

h4=plot(x,EE_PGA_TS_avr,altstyles{2});
% h5=plot(x,SW_apga.EE{2},altstyles{7});
% h6=plot(x,SW_npga.EE{1},altstyles{3});
% h7=plot(x,SW_npga.EE{2},altstyles{8});
h8=plot(x,EE_PGA_TbRS_avr,altstyles{4});
% h9=plot(x,SW_rdn.EE{2},altstyles{9});
h10=plot(x,EE_SW_rs_avr);


set(h10,'Color',[0,0.5,0],'Marker','d','LineStyle','-')%牛油绿
% set(h11,'Color',[0.6,0.2,0],'Marker','^','LineStyle','-')%brown
set(h0,'Color',[0.4940 0.1840 0.5560],'Marker','x','LineStyle','-')% purple


legend([h0,h1,h2,h3],'DBF','PS-HBF:4bits','PS-HBF:2bits','PS-HBF:1bit');
ylabel('Average EE (Gbits/J)')
xlabel('BSR')
xticks(x)
xlim([x(1) x(end)])
grid on 
box on

axesNew = axes('position',get(gca,'position'),'visible','off');
legend(axesNew,[h4,h8,h10],'SW-HBF:PGA-TS','SW-HBF:TS','SW-HBF:PGA-TbRs','location','west')


ccc=1;