function SE_sa=HBF_SA(params,Pb,Ns,N_rf,H_all,noise_power,episilon)
[Nr,Nt,Nc]=size(H_all);
Np=10;
power_avr=Pb/Ns;
F_ini=rand(Nt,N_rf);
F_ini_int = mapping(F_ini,Ns);
coef=power_avr/noise_power;
Ht=0;
for c=1:Nc
    Hc=H_all(:,:,c);
    HBF_tx.H_eff{c}=Hc'*Hc;
end

[SE_tx,F_rf]=SA_beamforming(HBF_tx,F_ini_int,coef,Ns,params);

Fp=F_rf'*F_rf;
[Uq,Sq]=eig(Fp);
Froot=Uq*sqrt(inv(Sq))*Uq';% (F'F)^(-0.5)
for c=1:Nc
    Hc=H_all(:,:,c);
    Qc=Hc*F_rf*Froot;
    [~,Sc,Vc]=svd(Qc);
    if N_rf>1
        rank_apx=min(Ns,trace(Sc(1:N_rf,:).^2>episilon));
    else
        rank_apx=Ns;
    end
    baseband_bf=Froot*Vc(:,1:rank_apx)*sqrt(power_avr);  
    precoder=F_rf*baseband_bf;
    HBF_info.transmitter_bf{c}=precoder;   
    HBF_info.H_eff{c}=Hc*precoder*(Hc*precoder)';
    
end
W_ini=rand(Nr,N_rf);
W_ini_int = mapping(W_ini,Ns);
coef=1/noise_power;
[SE_rx,W_rf]=SA_beamforming(HBF_info,W_ini_int,coef,Ns,params);
SE_sa=compute_rate(HBF_info,W_rf,coef);

end