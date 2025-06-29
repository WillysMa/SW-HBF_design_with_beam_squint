function [SE_ts,W_rf]=SW_HBF_random(Pb,Ns,N_rf,H_all,noise_power,episilon)


[Nr,Nt,Nc]=size(H_all);
Np=10;
power_avr=Pb/Ns;
F_ini=rand(Nt,N_rf);
F_rf = mapping(F_ini,Ns);

Fp=F_rf'*F_rf;
[Uq,Sq]=eig(Fp);
Froot=Uq*sqrt(inv(Sq))*Uq';% (F'F)^(-0.5)

for c=1:Nc
    Hc=H_all(:,:,c);

    Qc=Hc*F_rf*Froot;
    [~,Sc,Vc]=svd(Qc);
%     if N_rf>1
%         rank_apx=min(Ns,trace(Sc(1:N_rf,:).^2>episilon));
%     else
%         rank_apx=Ns;
%     end
    rank_apx=Ns;
    HBF_info.approx_rank{c}=rank_apx;
    baseband_bf=Froot*Vc(:,1:rank_apx)*sqrt(power_avr);  
    precoder=F_rf*baseband_bf;
    HBF_info.transmitter_bf{c}=precoder;   
    HBF_info.H_eff{c}=Hc*precoder*(Hc*precoder)';
    

end


W_ini=rand(Nr,N_rf);
W_rf = mapping(W_ini,Ns);
coef=1/noise_power;

SE_ts=compute_rate(HBF_info,W_rf,coef);


ccc=1;
end