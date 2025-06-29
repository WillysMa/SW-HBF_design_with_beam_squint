function SE=compute_rate_rx(HBF_info,BeaMtx,coef)
[~,Nc]=size(HBF_info.H_eff);
[~,N_rf]=size(BeaMtx);
rate=zeros(1,Nc);
for k=1:Nc
    rate(k)=real(log2(det(  eye(N_rf)+ coef* pinv(BeaMtx)*HBF_info.H_eff{k}*BeaMtx  )));
end
SE=mean(rate);
end