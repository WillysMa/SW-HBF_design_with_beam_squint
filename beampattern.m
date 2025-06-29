function beampattern(receive_bf,Ns,AtnDis,ofdm_paras)
Nc=ofdm_paras.Nc;% carrier number
D=ofdm_paras.D;% maximum delay
B=ofdm_paras.B;% bandwidth
fc=ofdm_paras.fc;% carrier frequency

% receive_bf=W_rf0_pg_ini1;
theta=0:0.01:2*pi;
len=length(theta);
[Nr,~]=size(receive_bf);
F_amp_all=zeros(Ns,len,Nc);
for c=1:Nc
%     f=fc;
    f=fc+(c-(Nc+1)/2)*B/Nc;
    array_response=steer_vector(theta,Nr,fc,f,AtnDis);
    for n=1:Ns
        arr_vec_ncns=receive_bf(:,n)'*array_response;
        F_amp_all(n,:,c)=abs(arr_vec_ncns)/sqrt(Nr);
       
    end
    figure
    pax = polaraxes;
    hold on
    for n=1:Ns
    polarplot(theta,F_amp_all(n,:,c),'LineWidth',1.2)
    end
%     legend('stream 1','stream 2','stream 3','stream 4')
%     title(['carrier=',num2str(c)])
    
    pax.ThetaDir = 'clockwise';
    pax.FontSize = 12;
     ccc=1;
end

end

%% defiend function
function array_vector=steer_vector(theta,N,fc,f,AtnDis)
len=length(theta);
a=zeros(N,len);
step=2*pi*AtnDis*sin(theta)*f/fc;
for n=1:N
    a(n,:)=exp(-1i*(n-1).*step);
end
array_vector=1/sqrt(N)*a;
end