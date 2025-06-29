
function [H,ch_pars]=Channel_Gen_compact(Nt,Nr,L,AtnDis,B,ofdm_paras)
Nc=ofdm_paras.Nc;% carrier number
D=ofdm_paras.D;% maximum delay
% B=ofdm_paras.B;% bandwidth
fc=ofdm_paras.fc;% carrier frequency
% L=10; % path number
Ts=1/B;
beta=1;
H=zeros(Nr,Nt,Nc);
var=1;
AoA_all=(2*rand(1,L)-1)*pi/2;ch_pars.AoA=AoA_all;
AoD_all=(2*rand(1,L)-1)*pi/2;ch_pars.AoD=AoD_all;
path_delay_all=rand(1,L)*(D-1)*Ts;ch_pars.path_delay=path_delay_all;
path_gain_all=sqrt(var/2)*(randn(L,Nc)+1i*randn(L,Nc));ch_pars.path_gain=path_gain_all;
Coef_matrix=zeros(L,Nc);
Path_Delay_gain=zeros(L,Nc);
for nc=1:Nc
    f=fc+(nc-(Nc+1)/2)*B/Nc;
    At=array_steering_dictionary(AoD_all,Nt,fc,f,AtnDis);
    Ar=array_steering_dictionary(AoA_all,Nr,fc,f,AtnDis);
    
    for l=1:L
        med=0;
        for d=1:D
         med=med+pulse_filter((d-1)*Ts-path_delay_all(l),Ts,beta)*exp(-1i*2*pi*(nc-1)*(d-1)/Nc);
        end
        Coef_matrix(l,nc)=med;
    end
    gain=sqrt(Nt*Nr/L)*Coef_matrix(:,nc).*path_gain_all(:,nc);
    H(:,:,nc)=Ar*diag(gain)*At';
    Path_Delay_gain(:,nc)=gain;
end       
ch_pars.Path_Delay_gain=Path_Delay_gain;
% ccc=1;
end

function p=pulse_filter(t,Ts,beta)
    if abs(t-Ts/2/beta)/abs(t)<1e-4 || abs(t+Ts/2/beta)/abs(t)<1e-4
        p=pi/4*sinc(1/2/beta);
    else
        p=sinc(t/Ts)*cos(pi*beta*t/Ts)/(1-(2*beta*t/Ts)^2);
    end
%     cc=1;
end

function array_matrix=array_steering_dictionary(angle_set,N,fc,f,AtnDis)
step=2*pi*AtnDis*sin(angle_set)*f/fc;
multiplier=0:N-1;
array_matrix=1/sqrt(N)*exp(-1i*multiplier'*step);

end
