function [SE_DBF,SE_DBFtx,DBF,Proj_Comb]=DBF_water_filling(Hc_all,P,Ns,noise_power)
[Nr,Nt,Nc]=size(Hc_all);
R_all=zeros(Nc,1);
rate_tx=zeros(Nc,1);
Proj_Comb=zeros(Nr,Nr,Nc);
Precoder_all=zeros(Nt,Ns,Nc);
Combiner_all=zeros(Nr,Ns,Nc);
for c=1:Nc
    Hc=Hc_all(:,:,c);
    [U,S,V]=svd(Hc);
sigular_vec=diag(S);
eigen_vec=sigular_vec.^2;

eigen_eff=eigen_vec(1:Ns);

power_allocatoin=water_filling(P,Ns,eigen_eff,noise_power);

power_matrix=diag(power_allocatoin(1:Ns));

U_effect=U(:,1:Ns);
V_effect=V(:,1:Ns);

precoder=V_effect*sqrt(power_matrix);
combiner=U_effect;

Precoder_all(:,:,c)=precoder;
Combiner_all(:,:,c)=combiner;
Proj_Comb(:,:,c)=combiner*pinv(combiner);

DBF.precoder{c}=precoder;
DBF.combiner{c}=combiner;


R_all(c)=real(log2(det( eye(Ns)+ 1/noise_power*pinv(combiner)*Hc*precoder*precoder'*Hc'*combiner  )));
rate_tx(c)=real(log2(det( eye(Nr)+ 1/noise_power*Hc*precoder*precoder'*Hc' )));
end
DBF.Precoder_all=Precoder_all;
DBF.Combiner_all=Combiner_all;
SE_DBF=mean(R_all);
SE_DBFtx=mean(rate_tx);
end
%% defined functions
function Rc=compute_rate(Hc,precoder,combiner,ch_rank,noise_power)
    noise_term=noise_power*(combiner'*combiner);
    other_term=combiner'*Hc*precoder*(combiner'*Hc*precoder)';
    det_term=eye(ch_rank)+inv(noise_term)*other_term;
    Rc=real(log2(det(det_term)));
end

