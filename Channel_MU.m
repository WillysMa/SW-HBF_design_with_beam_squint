function H=Channel_MU(Nt,Nu,L,AtnDis,B,ofdm_paras)

Nc=ofdm_paras.Nc;% carrier number
D=ofdm_paras.D;% maximum delay
% B=ofdm_paras.B;% bandwidth
fc=ofdm_paras.fc;% carrier frequency
% L=10; % path number
Ts=1/B;
beta=1;
H=zeros(Nt,Nu,Nc);
var=1;
AoD_all=(2*rand(Nu,L)-1)*pi/2;
path_delay_all=rand(Nu,L)*(D-1)*Ts;
path_gain_all=sqrt(var/2)*(randn(Nu,L)+1i*randn(Nu,L));
channel_power=0;
for nc=1:Nc
    f=fc+(nc-(Nc+1)/2)*B/Nc;
        for ii=1:Nu
            tmp_ch=0;
            for l=1:L
                path_gain=path_gain_all(ii,l);
                delay_gain=0;
                for d=1:D
                 delay_gain=delay_gain+pulse_filter((d-1)*Ts-path_delay_all(ii,l),Ts,beta)*exp(-1i*2*pi*(nc-1)*(d-1)/Nc);
                end
                array_vector=array_steering_dictionary(AoD_all(ii,l),Nt,fc,f,AtnDis);

                tmp_ch=tmp_ch+path_gain*delay_gain*array_vector;
            end 

             H(:,ii,nc)=sqrt(Nt/L)*tmp_ch;
            channel_power=channel_power+norm(H(:,ii,nc),'fro')^2;
        end   
end  
% disp(['Nt=',num2str(Nt),' channel power=', num2str(channel_power)])
% ccc=1;
end
%% defined functions
function p=pulse_filter(t,Ts,beta)
    if abs(t-Ts/2/beta)/abs(t)<1e-4 || abs(t+Ts/2/beta)/abs(t)<1e-4
        p=pi/4*sinc(1/2/beta);
    else
        p=sinc(t/Ts)*cos(pi*beta*t/Ts)/(1-(2*beta*t/Ts)^2);
    end
%     cc=1;
end

function array_vector=array_steering_dictionary(angle,N,fc,f,AtnDis)
step=2*pi*AtnDis*sin(angle)*f/fc;
multiplier=0:N-1;
array_vector=1/sqrt(N)*exp(-1i*multiplier'*step);

end
