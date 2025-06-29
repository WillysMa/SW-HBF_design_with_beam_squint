function [SE_HBF,W_rf,TC]=HBF_LSAA(Hc_all,P,Ns,N_rf,noise_power,episilon,phase_Set)
% 
tStart = tic;
[Nr,Nt,Nc]=size(Hc_all);
He=zeros(Nt,Nt);
for c=1:Nc
    Hc=Hc_all(:,:,c);
    He=He+Hc'* Hc;
end
H_avr=1/Nc*He;

power_avr=sqrt(P/(Ns));

coef_t=power_avr^2/noise_power;
V_rf=analog_bf_compute(coef_t,Nt,N_rf,H_avr);
V_rf=angle_qt(V_rf,phase_Set);

% check_max=max(abs(V_rf))
Q=V_rf'*V_rf;
[Uq,Sq]=eig(Q);
eigen_vec=diag(Sq);
eigen_s=1./sqrt(eigen_vec);
Sq_s=diag(eigen_s);
Q_s=Uq*Sq_s*Uq';
F_avr=zeros(Nr,Nr);

for c=1:Nc
    Hc=Hc_all(:,:,c);
    Hce=Hc*V_rf*Q_s;
%     rank_Hce=rank(Hce);
    [~,Sce,Vce]=svd(Hce);
    if N_rf>1
    eigen_ce=diag(Sce);
    rank_apx=min(Ns,sum(eigen_ce>episilon)); 
    else
        rank_apx=Ns;
    end
%     HBF_info.approx_rank{c}=rank_apx;
    baseband_bf=Q_s*Vce(:,1:rank_apx)*power_avr;
    precoder=V_rf*baseband_bf;
    HBF_info.transmitter_bf{c}=precoder;
    
    HBF_info.H_eff{c}=Hc*precoder*(Hc*precoder)';

    F_avr=F_avr+HBF_info.H_eff{c};
    
end
F_avr=1/Nc*F_avr;
%-----------------------------------combiner-------------------------------%
coef_r=1/noise_power/Nr;
W_rf=analog_bf_compute(coef_r,Nr,N_rf,F_avr);
W_rf=angle_qt(W_rf,phase_Set);

coef=1/noise_power;
SE_HBF=compute_rate(HBF_info,W_rf,coef);

TC=toc(tStart);
cc=1;
end
%% defined function



