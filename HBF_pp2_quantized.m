function [SE_HBF_pp,W_rf]=HBF_pp2_quantized(Hc_all,Pb,Ns,N_rf,noise_power,episilon,phase_Set)
[Nr,Nt,Nc]=size(Hc_all);
%-----------------------------precoder-------------------------------------%
Ht=0;
for c=1:Nc
    Hc=Hc_all(:,:,c);
    Rc=rank(Hc,1e-4);
    [~,~,V]=svd(Hc);
    Vc=V(:,1:Rc);
    Ht=Ht+Vc*Vc';
end
Ht=Ht/Nc;  

V_rf=analog_bf_compute(Ht,Nt,N_rf);
V_rf=angle_qt(V_rf,phase_Set);
% check_max=max(abs(V_rf))

Q=V_rf'*V_rf;
[Uq,Sq]=eig(Q);
eigen_vec=diag(Sq);
eigen_s=1./sqrt(eigen_vec);
Sq_s=diag(eigen_s);
Q_s=Uq*Sq_s*Uq';
Hr=0;
for c=1:Nc
    
    Hc=Hc_all(:,:,c);
    Hce=Hc*V_rf*Q_s;
%     rank_Hce=rank(Hce);
    [~,Sce,Vce]=svd(Hce);
    if N_rf>1
        sigular_ce=diag(Sce);
        rank_apx=min(Ns,sum(sigular_ce.^2>episilon));  
    else
        rank_apx=Ns;
    end
%     rank_apx=Ns;
    power_avr=sqrt(Pb/(rank_apx));
    
    HBF_info.approx_rank{c}=rank_apx;
    baseband_bf=Q_s*Vce(:,1:rank_apx)*power_avr;
    precoder=V_rf*baseband_bf;
    HBF_info.transmitter_bf{c}=precoder;
    HBF_info.H_eff{c}=Hc*precoder*(Hc*precoder)';
    
    Rc=rank(HBF_info.H_eff{c},1e-4);
    [~,~,V]=svd(HBF_info.H_eff{c});
    Vc=V(:,1:Rc);
    Hr=Hr+Vc*Vc';
    
end
Hr=Hr/Nc;
%-----------------------------------combiner-------------------------------%
W_rf=analog_bf_compute(Hr,Nr,N_rf);
W_rf=angle_qt(W_rf,phase_Set);
coef=1/noise_power;
SE_HBF_pp=compute_rate_rx(HBF_info,W_rf,coef);
ccc=1;
end

function analog_bf=analog_bf_compute(H_avr,N,N_rf)
% C=randn(N_rf,5*N_rf) + 1j*randn(N_rf,5*N_rf);
% B=C*C';% B can be any invertible matrix
% [U,~,~]=svd(H_avr);
Np=10;
U=block_power_method(H_avr,N_rf,Np);
Wb_ini=U;
W_rf=zeros(N,N_rf);
for row_id=1:N
    for col_id=1:N_rf
        W_rf(row_id,col_id)=exp(1j*angle(Wb_ini(row_id,col_id)));%1/sqrt(N)*
    end
end
analog_bf=W_rf;
end

