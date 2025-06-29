function H=channel_THz_MISO(Nt,Nu,AtnDis,B,ofdm_paras)
%  Note that the array gain is to be included for correctness.
radius=50;
Cell_CDT = zeros(6,2);%boundary of each cell
%% 正六边行小区内随机产生基站-用户位置
            Cell_CDT(1, :) = [-1, -1/sqrt(3)]*radius;
            Cell_CDT(2, :) = [-1,  1/sqrt(3)]*radius;
            Cell_CDT(3, :) = [ 0,  2/sqrt(3)]*radius;
            Cell_CDT(4, :) = [ 1,  1/sqrt(3)]*radius;
            Cell_CDT(5, :) = [ 1, -1/sqrt(3)]*radius;
            Cell_CDT(6, :) = [ 0, -2/sqrt(3)]*radius;
%基站位置选取
    BS_CDT = [0,0];
%     for j=1:Q
%         x_pos = BS_CDT(j, 1);
%         y_pos = BS_CDT(j, 2)*(2/sqrt(3) - (1/sqrt(3))*abs(x_pos));
%         BS_CDT(j,:) = [x_pos, y_pos]*radius;
%     end
%下行用户位置选取
    US_CDT = rand(Nu, 2)*2 - ones(Nu, 2);
    for j=1:Nu
        x_pos = US_CDT(j, 1);
        y_pos = US_CDT(j, 2)*(2/sqrt(3) - (1/sqrt(3))*abs(x_pos));
        US_CDT(j,:) = [x_pos, y_pos]*radius;
    end
%% 实际基站-用户距离计算 信道生成
%基站到下行用户信道
Dis_all=zeros(1,Nu);
for ii=1:Nu
    Dis_all(ii)=norm(BS_CDT-US_CDT(ii,:));
end

Nc=ofdm_paras.Nc;% carrier number
D=ofdm_paras.D;% maximum delay
% B=ofdm_paras.B;% bandwidth
fc=ofdm_paras.fc;% carrier frequency
Ts=1/B;
beta=1;
H=zeros(Nt,Nu,Nc);
MaCoef=0.0033;
light_speed=3*1e+8;

AoD_all=(2*rand(1,Nu)-1)*pi/2;
path_delay_all=rand(1,Nu)*(D-1)*Ts;

for nc=1:Nc
    f=fc+(nc-(Nc+1)/2)*B/Nc;
    for ii=1:Nu
        path_gain=light_speed/(4*pi*f*Dis_all(ii))*exp(-0.5*MaCoef*Dis_all(ii));% extremely small, needed to be compensated by array gain
        delay_gain=0;
        for d=1:D
         delay_gain=delay_gain+pulse_filter((d-1)*Ts-path_delay_all(ii),Ts,beta)*exp(-1i*2*pi*(nc-1)*(d-1)/Nc);
        end
        
        array_vector=array_steering_dictionary(AoD_all(ii),Nt,fc,f,AtnDis);
        H(:,ii,nc)=path_gain*delay_gain*array_vector;
    end
end       

ccc=1;
end

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
