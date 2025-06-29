clc;clear all;close all
G=1e+9;%bandwidth=1GHz
Nc=4;% carrier number
D=Nc/3;% maximum delay, see configure of simulation in 'Dynamic Subarrays for Hybrid Precoding inWideband mmWave MIMO Systems'
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



for nc=1:Nc
f=fc+(nc-(Nc+1)/2)*B/Nc;

% varepsilon=1+(nc/Nc-(Nc+1)/(2*Nc))*B/fc;

varepsilon=f/fc;
weight=randi(2,[Nt,1])-1;

theta=-1:0.001:1;
x=theta(1);

obj1=approx_array_gain(Nt,Nc,AtnDis,varepsilon,weight,x)

lenX=length(theta);
obj1=zeros(1,lenX);
obj2=zeros(1,lenX);
for ii=1:lenX

    
    obj1(ii)=array_gain_SW(weight,Nt,theta(ii),varepsilon);
    obj2(ii)=array_gain_SW(weight,Nt,-theta(ii),varepsilon);
end
dis=norm(obj2-obj1,2)
figure
subplot(121)
plot(obj1)

subplot(122)
plot(obj2)

ccc=1;



end
%%
function obj=approx_array_gain(Nt,Nc,AtnDis,varepsilon,weight,x)
Omega_tmp=zeros(Nt,Nc);
for c=1:Nc
    Omega_tmp(:,c)=exp(1i*2*AtnDis*pi*varepsilon*x*(1:Nt)');
end
Omega=Omega_tmp'/(Nc*sqrt(Nt));
obj=norm(Omega*weight,1)/sqrt(norm(weight,1));

end

function V=array_gain_SW(w,N,theta,varepsilon)
tmp=exp(1i*(1:N)*pi*varepsilon*theta);
V=abs(tmp*w);

end