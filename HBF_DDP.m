function [SE,SE1,TC]=HBF_DDP(Hc_all,Pb,Ns,N_rf,Ndt,Ndr,noise_power,AtnDis,B,ch_params,ofdm_paras,phase_Set)
tStart = tic;
episilon=1e-6;
[Nr,Nt,Nc]=size(Hc_all);
ofdm_paras.B=B;% bandwidth
% N_ttd=8;
[Au_all_tx, A_ttd_tx]=ttd_calculation(ch_params.AoD,Nt,Ndt, N_rf,Nc,AtnDis,ofdm_paras,phase_Set);
[Au_all_rx, A_ttd_rx]=ttd_calculation(ch_params.AoA,Nr,Ndr, N_rf,Nc,AtnDis,ofdm_paras,phase_Set);
rate1=zeros(1,Nc);   
coef=1/noise_power;
for nc=1:Nc
 
    Heff=Hc_all(:,:,nc)*Au_all_tx*A_ttd_tx(:,:,nc);
    [U,S,V]=svd(Heff);
    if N_rf>1
        sigular_ce=diag(S);
        rank_apx_tx=min(Ns,sum(sigular_ce>episilon));  
    else
        rank_apx_tx=Ns;
    end
%     rank_apx=Ns;
    power_avr=sqrt(Pb/(rank_apx_tx));
    digital_precoder=V(:,1:rank_apx_tx)*power_avr;
    HBF_info.digital_bf{nc}=digital_precoder;
    HBF_info.precoder{nc}=Au_all_tx*A_ttd_tx(:,:,nc)*digital_precoder;

    HBF_info.H_eff{nc}=Hc_all(:,:,nc)*HBF_info.precoder{nc}*(Hc_all(:,:,nc)*HBF_info.precoder{nc})';

    HBF_info.combiner{nc}=Au_all_rx*A_ttd_rx(:,:,nc);

   
%     rate1(nc)=real(log2(det(  eye(Nr)+ coef*HBF_info.H_eff{nc})));

    cccc=1;
                                       
end
SE1=mean(rate1);

rate2=zeros(1,Nc);
for k=1:Nc
    rate2(k)=real(log2(det(  eye(N_rf)+ coef* pinv(HBF_info.combiner{k})*HBF_info.H_eff{k}*HBF_info.combiner{k} )));
end
SE=mean(rate2);

TC=toc(tStart);
ccc=1;

end
%% defined functions
function [Au_all, A_ttd]=ttd_calculation(angle_set,N,N_ttd, N_rf,Nc,AtnDis,ofdm_paras,phase_Set)
fc=ofdm_paras.fc;% carrier frequency
B=ofdm_paras.B;
Ng=N/N_ttd;
Au_all=[];
Delay=zeros(N_ttd,N_rf);
for ll=1:N_rf
    angle_eff=sin(angle_set(ll));
    xx=Ng*angle_eff/2;            
    array_vector=array_steering_dictionary(angle_eff,N,fc,fc,AtnDis);    
    
    Au=zeros(N,N_ttd);
    for kk=1:N_ttd
        if angle_eff<0
            Delay(kk,ll)=(N_ttd-1)*abs(xx)/fc+kk*xx/fc;
        else
            Delay(kk,ll)=kk*xx/fc;
        end
        Au((kk-1)*Ng+1:kk*Ng,kk)=array_vector((kk-1)*Ng+1:kk*Ng)*exp(1i*pi*(kk-1)*Ng*angle_eff);
    end
    Au_all=[Au_all Au];
end    
Au_all=angle_qt(Au_all,phase_Set);% low-resolution PSs

A_ttd=zeros(N_rf*N_ttd,N_rf,Nc);
for nc=1:Nc
    f=fc+(nc-(Nc+1)/2)*B/Nc;
    for ll=1:N_rf
        A_ttd((ll-1)*N_ttd+1:ll*N_ttd,ll,nc)=exp(-1i*2*pi*f*Delay(:,ll));
    end
end


end

function array_matrix=array_steering_dictionary(angle_eff,N,fc,f,AtnDis)
step=2*pi*AtnDis*angle_eff*f/fc;
multiplier=0:N-1;
array_matrix=1/sqrt(N)*exp(-1i*multiplier'*step);


end