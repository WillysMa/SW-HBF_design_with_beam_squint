function SE=compute_rate_v0(data,BeaMtx,coef)
[~,Nc]=size(data.H_eff);
[U,~]=qr(BeaMtx,0);
rate=zeros(1,Nc);

for k=1:Nc
    Surrgate=U'*data.H_eff{k}*U;
    Eigen_all=eig(Surrgate);   
    rate(k)=real(sum(log2(  1+ coef* Eigen_all  )));
end
SE=mean(rate);
end