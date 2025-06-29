clc;clear all;close all
GHz=1e+9;%bandwidth=1GHz
Nc=60;% carrier number
D=Nc/3;% maximum delay, see configure of simulation in 'Dynamic Subarrays for Hybrid Precoding inWideband mmWave MIMO Systems'
B=30*GHz;% bandwidth
fc=300*GHz;% carrier frequency
AtnDis=1/2;
Nt=16;
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

BF_sw=ones(Nt,1);

Ball=(5:5:30)*GHz;
lenB=length(Ball);
Score_ps_true=zeros(4,lenB);
Score_ps_apx=zeros(4,lenB);


% weight=[1 1 0 1 1 1 0 0 1 1 0 1 1 0 1 1]';

% theta_all=-1:0.001:1;
% lenN=length(theta_all);
% output=zeros(1,lenN);
% for ii=1:lenN
% output(ii)=obj_SW_array(ofdm_paras,AtnDis,weight,theta_all(ii));
% end
% figure
% plot(theta_all,output)

Score_sw_true=zeros(4,lenB);
Score_sw_apx=zeros(4,lenB);
BSR_all=zeros(1,lenB);
for i=1:lenB
    ofdm_paras.B=Ball(i);
    Nt=16;
    [Score_ps_true(1,i),Score_ps_apx(1,i),BSR_all(i)]=PS_array_gain(ofdm_paras,Nt,AtnDis);
    weight=randi(2,[Nt,1])-1;
    [Score_sw_true(1,i),Score_sw_apx(1,i)]=SW_array_gain(ofdm_paras,AtnDis,weight,Nt);
    Nt=64;
    weight=randi(2,[Nt,1])-1;
    [Score_ps_true(2,i),Score_ps_apx(2,i),BSR_all(i)]=PS_array_gain(ofdm_paras,Nt,AtnDis);
    [Score_sw_true(2,i),Score_sw_apx(2,i)]=SW_array_gain(ofdm_paras,AtnDis,weight,Nt);
    Nt=256;
    weight=randi(2,[Nt,1])-1;
    [Score_ps_true(3,i),Score_ps_apx(3,i),BSR_all(i)]=PS_array_gain(ofdm_paras,Nt,AtnDis);
    [Score_sw_true(3,i),Score_sw_apx(3,i)]=SW_array_gain(ofdm_paras,AtnDis,weight,Nt);
    
    Nt=4096;
    weight=randi(2,[Nt,1])-1;
    [Score_ps_true(4,i),Score_ps_apx(4,i),BSR_all(i)]=PS_array_gain(ofdm_paras,Nt,AtnDis);
    [Score_sw_true(4,i),Score_sw_apx(4,i)]=SW_array_gain(ofdm_paras,AtnDis,weight,Nt);

%     Score_sw(i)=2/3;
end
figure
x=Ball/GHz;
hold on
h1=plot(x,Score_ps_true(1,:),'r-s',x,Score_ps_true(2,:),'b-o',x,Score_ps_true(3,:),'k-p',...
    x,Score_ps_true(4,:),'m-<','LineWidth',2)
h2=plot(x,Score_ps_apx(1,:),'r:s',x,Score_ps_apx(2,:),'b:o',x,Score_ps_apx(3,:),'k:p',...
    x,Score_ps_apx(4,:),'m:<','LineWidth',2)
legend([h1(1),h1(2),h1(3),h1(4),h2(1),h2(2),h2(3),h2(4)],'True $N=16$','True $N=64$','True $N=256$','True $N=4096$', ...
    'Apx $N=16$','Apx $N=64$','Apx $N=256$','Apx $N=4096$','interpreter','latex')
xlabel('Bandwidth(GHz)')
ylabel('EAG-PS')
grid on

figure
x=Ball/GHz;
hold on
h1=plot(x,Score_sw_true(1,:),'r-s',x,Score_sw_true(2,:),'b-o',x,Score_sw_true(3,:),'k-p',...
    x,Score_ps_true(4,:),'m-<','LineWidth',2)
h2=plot(x,Score_sw_apx(1,:),'r:s',x,Score_sw_apx(2,:),'b:o',x,Score_sw_apx(3,:),'k:p',...
    x,Score_sw_apx(4,:),'m:<','LineWidth',2)
legend([h1(1),h1(2),h1(3),h1(4),h2(1),h2(2),h2(3),h2(4)],'True $N=16$','True $N=64$','True $N=256$','True $N=4096$', ...
    'Apx $N=16$','Apx $N=64$','Apx $N=256$','Apx $N=4096$','interpreter','latex')
xlabel('Bandwidth(GHz)')
ylabel('EAG-SW')
grid on
ccc=1;

%% defined function
function obj=obj_SW_array(ofdm_paras,AtnDis,weight,theta)
[Nt,~]=size(weight);
B=ofdm_paras.B;
Nc=ofdm_paras.Nc;
fc=ofdm_paras.fc;
Omega_tmp=zeros(Nt,Nc);
for nc=1:Nc
    f=fc+(nc-(Nc+1)/2)*B/Nc;
    varepsilon=f/fc;
    Omega_tmp(:,nc)=exp(1i*2*AtnDis*pi*varepsilon*theta*(1:Nt)');
end
Omega=Omega_tmp'/(Nc*sqrt(Nt));
obj=norm(Omega*weight,1)/sqrt(norm(weight,1));
end

function [Score_ture,Score_apx]=SW_array_gain(ofdm_paras,AtnDis,weight,Nt)
M=10000;
temp=0;
for m=1:M
    theta=m/M;
    obj=obj_SW_array(ofdm_paras,AtnDis,weight,theta);
    temp=temp+obj;
end
Score_ture=temp/M;
Score_apx=1/3*obj_SW_array(ofdm_paras,AtnDis,weight,1)+2/3*sqrt(norm(weight,1)/Nt);

end